use super::{
    cache::{Cache, Object},
    client::download_url,
    *,
};
use md5::Md5;
use reqwest::{
    blocking::Client,
    header::{CONTENT_RANGE, ETAG, IF_RANGE, RANGE},
    StatusCode,
};
use sha1::Sha1;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::Path,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    },
    time::{Duration, Instant},
};

pub const DEFAULT_SAFETY: u64 = 1_000_000_000;
#[derive(Debug, Clone)]
pub struct Budgets {
    pub max_download: u64,
    pub max_working_space: Option<u64>,
    pub safety: u64,
}
impl Default for Budgets {
    fn default() -> Self {
        Self {
            max_download: 1_000_000_000,
            max_working_space: None,
            safety: DEFAULT_SAFETY,
        }
    }
}
#[derive(Debug, serde::Serialize)]
pub struct Plan {
    pub accession: Pxd,
    pub selected_files: Vec<(String, String, u64)>,
    pub download_bytes: u64,
    pub cache_current_bytes: u64,
    pub cache_limit_bytes: u64,
    pub free_filesystem_bytes: u64,
    pub temporary_workspace_bytes: u64,
    pub derived_artifact_bytes: u64,
    pub result_budget_bytes: u64,
    pub safety_margin_bytes: u64,
    pub peak_working_bytes: u64,
    pub expected_evictions: Vec<String>,
    pub expected_final_large_data_bytes: u64,
    pub expected_retained_results_provenance: String,
}
pub fn plan(
    cache: &Cache,
    m: &Manifest,
    files: &[RemoteFile],
    budgets: &Budgets,
    ephemeral: bool,
    batch_size: usize,
    result_budget: u64,
) -> Result<Plan> {
    if files.is_empty() {
        return Err("no files selected; inspect `pride files` and choose --file, --format, --category, or --all".into());
    }
    if batch_size == 0 {
        return Err("batch size must be positive".into());
    }
    let mut unique = BTreeMap::new();
    for f in files {
        download_url(f)?;
        unique.entry(f.object_key()).or_insert(f.size()?);
    }
    let status = cache.status()?;
    let mut missing = BTreeMap::new();
    let mut additional = BTreeMap::new();
    for f in files {
        let key = f.object_key();
        let size = f.size()?;
        let (present, reusable, resumable) = cached_bytes(cache, f)?;
        missing.insert(
            key.clone(),
            size.saturating_sub(if reusable { present } else { resumable }),
        );
        additional.insert(key, size.saturating_sub(present));
    }
    let download = total(missing.values().copied())?;
    if download > budgets.max_download {
        return Err(format!(
            "download budget exceeded: {download} bytes required, --max-download {}",
            budgets.max_download
        )
        .into());
    }
    let release = ephemeral
        && !cache.index.pinned.contains(&m.accession)
        && !unique.keys().any(|key| {
            cache
                .index
                .objects
                .get(key)
                .is_some_and(|o| cache.pinned(o) || o.retention == Retention::Keep)
        });
    let batch_peak = |values: Vec<u64>| -> Result<u64> {
        let mut sizes = values;
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        total(
            sizes
                .into_iter()
                .take(if release { batch_size } else { usize::MAX }),
        )
    };
    let working_sources = batch_peak(unique.values().copied().collect())?;
    let held = batch_peak(additional.values().copied().collect())?;
    let result_workspace = result_budget / files.len() as u64;
    let managed_peak = total([working_sources, result_workspace])?;
    let peak = total([working_sources, result_budget, budgets.safety])?;
    if budgets.max_working_space.is_some_and(|b| peak > b) {
        return Err(format!("working-space budget exceeded: need {peak} bytes including output allowance and safety margin").into());
    }
    if managed_peak > cache.limit {
        return Err(format!("selected working set including temporary results needs {managed_peak} bytes, exceeding hard cache ceiling {}; reduce --batch-size or select smaller files",cache.limit).into());
    }
    let exclude = unique.keys().cloned().collect();
    let evictions = cache.eviction_plan(
        total([held, result_workspace])?,
        total([held, result_budget, budgets.safety])?,
        &exclude,
    )?;
    let freed = total(
        evictions
            .iter()
            .map(|k| cache.index.objects.get(k).unwrap())
            .map(|o| {
                cache
                    .path(&o.relative_path)
                    .and_then(|p| Ok(fs::metadata(p)?.len()))
            })
            .collect::<Result<Vec<_>>>()?,
    )?;
    let retained = if release {
        let selected_present = total(
            unique
                .keys()
                .filter_map(|key| cache.index.objects.get(key))
                .filter_map(|o| cache.path(&o.relative_path).ok())
                .filter_map(|p| fs::metadata(p).ok())
                .map(|m| m.len()),
        )?;
        status
            .large_data_bytes
            .saturating_sub(freed)
            .saturating_sub(selected_present)
    } else {
        total([status.large_data_bytes, held])?.saturating_sub(freed)
    };
    Ok(Plan { accession:m.accession.clone(),selected_files:files.iter().map(|f|Ok((f.id.clone(),f.filename.clone(),f.size()?))).collect::<Result<_>>()?,
        download_bytes:download,cache_current_bytes:status.large_data_bytes,cache_limit_bytes:cache.limit,free_filesystem_bytes:status.free_filesystem_bytes,
        temporary_workspace_bytes:result_workspace,derived_artifact_bytes:0,result_budget_bytes:result_budget,safety_margin_bytes:budgets.safety,peak_working_bytes:peak,
        expected_evictions:evictions,expected_final_large_data_bytes:retained,expected_retained_results_provenance:"Metadata, source identities, checksums, configurations and verified final results; exact output size is unknown before analysis".into() })
}

// Separate bytes physically present from bytes safely reusable/resumable. Partial
// data already occupies the cache and must not be reserved a second time.
fn cached_bytes(cache: &Cache, file: &RemoteFile) -> Result<(u64, bool, u64)> {
    let Some(o) = cache.index.objects.get(&file.object_key()) else {
        return Ok((0, false, 0));
    };
    let path = cache.path(&o.relative_path)?;
    if !path.is_file() {
        return Ok((0, false, 0));
    }
    let size = file.size()?;
    let present = fs::metadata(path)?.len();
    let reusable = matches!(o.state, State::Verified | State::DownloadedUnverified)
        && o.local_sha256.is_some()
        && present == size;
    let resumable = if matches!(o.state, State::Partial | State::Downloading)
        && o.relative_path.starts_with("tmp/")
        && present <= size
        && (!file.checksums.is_empty() || o.etag.is_some())
    {
        present
    } else {
        0
    };
    Ok((present, reusable, resumable))
}

pub struct Downloader {
    client: Client,
    pub cancelled: Arc<AtomicBool>,
    pub remaining_download: u64,
    pub safety: u64,
}
impl Downloader {
    pub fn new(budgets: &Budgets, cancelled: Arc<AtomicBool>) -> Result<Self> {
        Ok(Self {
            client: Client::builder()
                .connect_timeout(Duration::from_secs(15))
                .timeout(Duration::from_secs(3600))
                .user_agent(concat!("percolator-rs/", env!("CARGO_PKG_VERSION")))
                .redirect(reqwest::redirect::Policy::custom(|a| {
                    if a.previous().len() >= 5
                        || (a.url().scheme() != "https" && a.url().host_str() != Some("127.0.0.1"))
                    {
                        a.stop()
                    } else {
                        a.follow()
                    }
                }))
                .build()?,
            cancelled,
            remaining_download: budgets.max_download,
            safety: budgets.safety,
        })
    }
    pub fn fetch(
        &mut self,
        cache: &mut Cache,
        m: &mut Manifest,
        f: &RemoteFile,
        protected: &BTreeSet<String>,
    ) -> Result<String> {
        cache.require_write()?;
        let size = f.size()?;
        let key = f.object_key();
        if f.untyped_checksum.is_some() {
            return Err(format!(
                "{} has an authoritative checksum with unknown algorithm; refusing to bypass it",
                f.filename
            )
            .into());
        }
        if let Some(existing) = cache.index.objects.get(&key).cloned() {
            let path = cache.path(&existing.relative_path)?;
            if matches!(
                existing.state,
                State::Verified | State::DownloadedUnverified | State::Prepared
            ) && path.is_file()
            {
                let hashes = hash_file(&path)?;
                let verification = verify(f, &hashes);
                if hashes.bytes != size
                    || existing.local_sha256.as_deref() != Some(hashes.sha256.as_str())
                    || verification.as_ref().is_err()
                {
                    let mut o = existing;
                    o.state = State::Corrupt;
                    o.verification = verification_records(f, &hashes);
                    cache.record_download(m, f, o)?;
                    return Err(format!(
                        "cached {} is corrupt; validity revoked. Retry fetch to replace it",
                        f.filename
                    )
                    .into());
                }
                let mut o = existing;
                o.projects.insert(m.accession.clone());
                o.last_used_unix_seconds = now();
                o.verification = verification?;
                cache.record_download(m, f, o)?;
                eprintln!("reuse {:?}: local SHA-256 checked", f.filename);
                return Ok(key);
            }
        }
        let (_, _, resumable) = cached_bytes(cache, f)?;
        let needed_transfer = size.saturating_sub(resumable);
        if self.remaining_download < needed_transfer {
            return Err(format!("remaining transfer budget {} cannot cover {} ({needed_transfer} bytes remaining); retries/restarts also consume the operation budget",self.remaining_download,f.filename).into());
        }
        // Reserve the entire object before receiving data, including any existing partial.
        let extra = cache
            .index
            .objects
            .get(&key)
            .map(|o| cache.path(&o.relative_path))
            .transpose()?
            .and_then(|p| fs::metadata(p).ok())
            .map(|x| x.len())
            .unwrap_or(0);
        let evictions = cache.eviction_plan(
            size.saturating_sub(extra),
            total([size.saturating_sub(extra), self.safety])?,
            protected,
        )?;
        cache.evict(&evictions, false)?;
        let part_rel = format!("tmp/{key}.part");
        let part = cache.path(&part_rel)?;
        let mut o = cache.index.objects.get(&key).cloned().unwrap_or(Object {
            key: key.clone(),
            relative_path: part_rel.clone(),
            bytes: size,
            state: State::Remote,
            local_sha256: None,
            verification: vec![],
            projects: BTreeSet::new(),
            last_used_unix_seconds: now(),
            retention: Retention::KeepIfPinned,
            result_verified: false,
            reproducible: true,
            etag: None,
            pin_validated: false,
        });
        if o.relative_path != part_rel {
            let old = cache.path(&o.relative_path)?;
            if old.is_file() {
                fs::remove_file(old)?;
            }
        }
        if o.state == State::Corrupt && part.exists() {
            fs::remove_file(&part)?;
        }
        o.relative_path = part_rel;
        o.state = State::Downloading;
        o.pin_validated = false;
        o.local_sha256 = None;
        o.verification.clear();
        o.bytes = size;
        o.projects.insert(m.accession.clone());
        cache.record_download(m, f, o.clone())?;
        let outcome = self.transfer(cache, m, f, &mut o, &part);
        match outcome {
            Ok(hashes) => {
                o.local_sha256 = Some(hashes.sha256.clone());
                o.verification = verification_records(f, &hashes);
                if hashes.bytes != size || o.verification.iter().any(|v| !v.verified) {
                    o.state = State::Corrupt;
                    cache.record_download(m, f, o)?;
                    return Err(format!("checksum/size mismatch for {:?}; corrupt partial preserved, never published",f.filename).into());
                }
                o.state = if o.verification.is_empty() {
                    State::DownloadedUnverified
                } else {
                    State::Verified
                };
                let final_rel = format!("objects/{key}");
                let dest = cache.path(&final_rel)?;
                if dest.exists() {
                    return Err("unexpected destination collision; refusing overwrite".into());
                }
                fs::rename(&part, &dest)?;
                File::open(dest.parent().unwrap())?.sync_all()?;
                o.relative_path = final_rel;
                o.last_used_unix_seconds = now();
                cache.record_download(m, f, o)?;
                eprintln!(
                    "completed {:?}: {} bytes; {}",
                    f.filename,
                    size,
                    if f.checksums.is_empty() {
                        "repository verification unavailable; local SHA-256 recorded"
                    } else {
                        "repository checksums VERIFIED; local SHA-256 recorded"
                    }
                );
                Ok(key)
            }
            Err(e) => {
                o.state = State::Partial;
                cache.record_download(m, f, o)?;
                Err(e)
            }
        }
    }
    fn transfer(
        &mut self,
        cache: &mut Cache,
        m: &mut Manifest,
        f: &RemoteFile,
        o: &mut Object,
        part: &Path,
    ) -> Result<Hashes> {
        let url = download_url(f)?;
        let size = f.size()?;
        let mut last = String::new();
        for attempt in 0..3 {
            if self.cancelled.load(Ordering::Relaxed) {
                return Err("download interrupted; partial is resumable or safely prunable".into());
            }
            let mut offset = fs::metadata(part).map(|x| x.len()).unwrap_or(0);
            if offset > size || (offset > 0 && f.checksums.is_empty() && o.etag.is_none()) {
                File::create(part)?.sync_all()?;
                offset = 0;
            }
            if offset == size && part.is_file() {
                return hash_file(part);
            }
            let mut request = self
                .client
                .get(url.clone())
                .header("Accept-Encoding", "identity");
            if offset > 0 {
                request = request.header(RANGE, format!("bytes={offset}-"));
                if let Some(tag) = &o.etag {
                    request = request.header(IF_RANGE, tag);
                }
            }
            let response = request.send();
            let mut response = match response {
                Ok(r) if r.status().is_success() => r,
                Ok(r) => {
                    last = format!("download HTTP {} for {:?}", r.status(), f.filename);
                    if !r.status().is_server_error() && r.status().as_u16() != 429 {
                        return Err(last.into());
                    }
                    if attempt < 2 {
                        std::thread::sleep(Duration::from_millis(250 << attempt));
                    }
                    continue;
                }
                Err(e) => {
                    last = e.to_string();
                    if attempt < 2 {
                        std::thread::sleep(Duration::from_millis(250 << attempt));
                    }
                    continue;
                }
            };
            if response.status() == StatusCode::PARTIAL_CONTENT {
                let expected = format!("bytes {}-{}/{}", offset, size.saturating_sub(1), size);
                if response
                    .headers()
                    .get(CONTENT_RANGE)
                    .and_then(|v| v.to_str().ok())
                    != Some(expected.as_str())
                {
                    return Err("invalid Content-Range; refusing unsafe resume".into());
                }
            } else if response.status() == StatusCode::OK {
                offset = 0;
            } else {
                return Err("unexpected successful download HTTP status".into());
            }
            if response
                .content_length()
                .is_some_and(|n| n != size - offset)
            {
                return Err("remote Content-Length disagrees with PRIDE size; refresh metadata before retrying".into());
            }
            if let Some(tag) = response
                .headers()
                .get(ETAG)
                .and_then(|x| x.to_str().ok())
                .filter(|x| !x.starts_with("W/"))
            {
                o.etag = Some(tag.into());
            }
            cache.record_download(m, f, o.clone())?;
            let mut out = OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(offset == 0)
                .append(offset > 0)
                .open(part)?;
            let mut buf = [0u8; 64 * 1024];
            let mut progress = Instant::now();
            let mut pos = offset;
            let transfer = (|| -> Result<()> {
                loop {
                    if self.cancelled.load(Ordering::Relaxed) {
                        return Err("download interrupted".into());
                    }
                    if pos == size {
                        let mut extra = [0u8; 1];
                        if response.read(&mut extra)? != 0 {
                            return Err("remote response exceeds declared size".into());
                        }
                        break;
                    }
                    if self.remaining_download == 0 {
                        return Err("download transfer budget exhausted".into());
                    }
                    let nmax = (buf.len() as u64)
                        .min(size - pos)
                        .min(self.remaining_download) as usize;
                    let n = response.read(&mut buf[..nmax])?;
                    if n == 0 {
                        return Err("interrupted/truncated download body".into());
                    }
                    self.remaining_download -= n as u64;
                    // Check actual disk space each chunk; managed capacity was reserved
                    // under the global lock and the response is bounded by declared size.
                    if cache::available_space(&cache.root)? < total([n as u64, self.safety])? {
                        return Err("filesystem free-space safety margin reached".into());
                    }
                    out.write_all(&buf[..n])?;
                    pos += n as u64;
                    if progress.elapsed() > Duration::from_secs(2) {
                        eprintln!("download {:?}: {pos}/{size} bytes", f.filename);
                        progress = Instant::now();
                    }
                }
                out.sync_all()?;
                Ok(())
            })();
            match transfer {
                Ok(()) => return hash_file(part),
                Err(e) => {
                    out.sync_all()?;
                    last = e.to_string();
                    if self.cancelled.load(Ordering::Relaxed)
                        || last.contains("budget")
                        || last.contains("space")
                    {
                        return Err(e);
                    }
                }
            }
        }
        Err(format!("download failed after three attempts: {last}; partial retained").into())
    }
}
#[derive(Debug)]
pub struct Hashes {
    pub bytes: u64,
    pub md5: String,
    pub sha1: String,
    pub sha256: String,
}
pub fn hash_file(path: &Path) -> Result<Hashes> {
    let mut file = File::open(path)?;
    let mut m = Md5::new();
    let mut s = Sha1::new();
    let mut h = Sha256::new();
    let mut bytes = 0;
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        m.update(&buf[..n]);
        s.update(&buf[..n]);
        h.update(&buf[..n]);
        bytes += n as u64;
    }
    Ok(Hashes {
        bytes,
        md5: format!("{:x}", m.finalize()),
        sha1: format!("{:x}", s.finalize()),
        sha256: format!("{:x}", h.finalize()),
    })
}
pub fn verification_records(file: &RemoteFile, hashes: &Hashes) -> Vec<Verification> {
    file.checksums
        .iter()
        .map(|c| {
            let actual = match c.algorithm {
                Algorithm::Md5 => &hashes.md5,
                Algorithm::Sha1 => &hashes.sha1,
                Algorithm::Sha256 => &hashes.sha256,
            };
            Verification {
                algorithm: c.algorithm,
                expected: c.value.clone(),
                actual: actual.clone(),
                verified: c.value == *actual,
                authority: c.authority.clone(),
            }
        })
        .collect()
}
pub fn verify(file: &RemoteFile, hashes: &Hashes) -> Result<Vec<Verification>> {
    let v = verification_records(file, hashes);
    if v.iter().any(|v| !v.verified) {
        Err("authoritative checksum mismatch".into())
    } else {
        Ok(v)
    }
}
