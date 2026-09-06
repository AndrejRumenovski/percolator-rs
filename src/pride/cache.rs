//! Serialized, bounded, disposable data cache. Metadata/results are never pruning targets.
use super::*;
use fs2::FileExt;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Component, Path, PathBuf},
};

pub const DEFAULT_LIMIT: u64 = 50_000_000_000;
const MARKER: &str = ".percolator-pride-cache-v1";
const OWNER: &str = "percolator-rs PRIDE cache schema 1\n";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Object {
    pub key: String,
    pub relative_path: String,
    pub bytes: u64,
    pub state: State,
    pub local_sha256: Option<String>,
    pub verification: Vec<Verification>,
    pub projects: BTreeSet<Pxd>,
    pub last_used_unix_seconds: u64,
    pub retention: Retention,
    pub result_verified: bool,
    pub reproducible: bool,
    pub etag: Option<String>,
    pub pin_validated: bool,
}
impl Object {
    pub fn local(&self) -> LocalFile {
        LocalFile {
            object_key: self.key.clone(),
            local_relative_path: self.relative_path.clone(),
            state: self.state,
            local_sha256: self.local_sha256.clone(),
            verification: self.verification.clone(),
            verification_unavailable: self.verification.is_empty(),
            pin_validated: self.pin_validated,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Index {
    pub schema_version: u32,
    pub pinned: BTreeSet<Pxd>,
    pub objects: BTreeMap<String, Object>,
}
impl Default for Index {
    fn default() -> Self {
        Self {
            schema_version: 1,
            pinned: BTreeSet::new(),
            objects: BTreeMap::new(),
        }
    }
}
#[derive(Debug, Serialize)]
pub struct Status {
    pub root: PathBuf,
    pub limit_bytes: u64,
    pub large_data_bytes: u64,
    pub source_bytes: u64,
    pub prepared_bytes: u64,
    pub temporary_partial_bytes: u64,
    pub pinned_bytes: u64,
    pub evictable_bytes: u64,
    pub protected_by_retention_bytes: u64,
    pub untracked_bytes: u64,
    pub metadata_bytes: u64,
    pub results_bytes: u64,
    pub project_count: usize,
    pub free_filesystem_bytes: u64,
    pub largest_objects: Vec<(String, u64)>,
}
#[derive(Debug, Serialize)]
pub struct Cleanup {
    pub before_bytes: u64,
    pub objects: Vec<String>,
    pub freed_bytes: u64,
    pub remaining_bytes: u64,
    pub pinned_remaining_bytes: u64,
    pub provenance_results_retained: bool,
    pub dry_run: bool,
}

pub struct Cache {
    pub root: PathBuf,
    pub limit: u64,
    pub index: Index,
    writable: bool,
    _lock: Option<File>,
}
impl Drop for Cache {
    fn drop(&mut self) {
        // Explicit unlock also releases a shared open-file description briefly
        // inherited by a concurrently forked child before its close-on-exec runs.
        // Closing only this descriptor can otherwise leave a phantom lock behind.
        if let Some(lock) = &self._lock {
            let _ = FileExt::unlock(lock);
        }
    }
}
impl Cache {
    pub fn default_root() -> Result<PathBuf> {
        if let Some(p) = std::env::var_os("PERCOLATOR_RS_PRIDE_CACHE") {
            return Ok(p.into());
        }
        if let Some(p) = std::env::var_os("XDG_CACHE_HOME") {
            return Ok(PathBuf::from(p).join("percolator-rs/pride"));
        }
        Ok(
            PathBuf::from(std::env::var_os("HOME").ok_or("set PERCOLATOR_RS_PRIDE_CACHE or HOME")?)
                .join(".cache/percolator-rs/pride"),
        )
    }
    /// Read-only mode creates nothing, including when the requested cache does not exist.
    pub fn open(root: &Path, limit: u64, writable: bool) -> Result<Self> {
        if limit == 0 {
            return Err("cache ceiling must be greater than zero".into());
        }
        let root = if root.is_absolute() {
            root.to_path_buf()
        } else {
            std::env::current_dir()?.join(root)
        };
        reject_symlinks(&root)?;
        for p in root.ancestors() {
            if p.join(".git").exists() {
                return Err("PRIDE large-data cache must be outside a Git working tree; set --cache-dir or PERCOLATOR_RS_PRIDE_CACHE".into());
            }
        }
        if !root.exists() && !writable {
            return Ok(Self {
                root,
                limit,
                index: Index::default(),
                writable,
                _lock: None,
            });
        }
        if !root.join(MARKER).exists() {
            if root.exists() && fs::read_dir(&root)?.next().is_some() {
                return Err(
                    "refusing to claim a nonempty directory without PRIDE ownership marker".into(),
                );
            }
            if !writable {
                return Err("directory is not an initialized PRIDE cache".into());
            }
            fs::create_dir_all(&root)?;
            let mut f = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(root.join(MARKER))?;
            f.write_all(OWNER.as_bytes())?;
            f.sync_all()?;
        }
        reject_symlinks(&root.join(MARKER))?;
        if fs::read_to_string(root.join(MARKER))? != OWNER {
            return Err("unrecognized PRIDE ownership marker".into());
        }
        let lockpath = root.join("cache.lock");
        reject_symlinks(&lockpath)?;
        let lock = if writable {
            let f = OpenOptions::new()
                .create(true)
                .truncate(false)
                .read(true)
                .write(true)
                .open(lockpath)?;
            f.try_lock_exclusive().map_err(|_| {
                "PRIDE cache is in use by another operation; try again after it finishes"
            })?;
            Some(f)
        } else if lockpath.exists() {
            let f = File::open(lockpath)?;
            FileExt::try_lock_shared(&f)
                .map_err(|_| "PRIDE cache is in use by another operation")?;
            Some(f)
        } else {
            None
        };
        for dir in ["objects", "prepared", "tmp", "manifests", "results"] {
            let p = root.join(dir);
            reject_symlinks(&p)?;
            if writable {
                fs::create_dir_all(p)?;
            }
        }
        let index_path = root.join("index.json");
        reject_symlinks(&index_path)?;
        let index = if index_path.exists() {
            read_json(&index_path).map_err(|e| format!("corrupt PRIDE index; no data deleted: {e}. Restore index.json from backup; manifests and results remain intact."))?
        } else {
            if ["objects", "prepared", "tmp"].iter().any(|d| {
                fs::read_dir(root.join(d))
                    .ok()
                    .is_some_and(|mut r| r.next().is_some())
            }) {
                return Err("missing PRIDE index with existing data; refusing to infer validity or discard pinning".into());
            }
            Index::default()
        };
        let mut c = Self {
            root,
            limit,
            index,
            writable,
            _lock: lock,
        };
        c.validate_index()?;
        // Reject nested symlinks/special files before any destructive operation.
        c.status()?;
        if writable {
            c.recover_interrupted()?;
        }
        Ok(c)
    }
    fn recover_interrupted(&mut self) -> Result<()> {
        let mut changed = false;
        for key in self.index.objects.keys().cloned().collect::<Vec<_>>() {
            let object = &self.index.objects[&key];
            let recorded = self.path(&object.relative_path)?;
            let folder = if key.starts_with("prepared-") {
                "prepared"
            } else {
                "objects"
            };
            let completed = self.path(&format!("{folder}/{key}"))?;
            // Rename can succeed before the index commit. Recover ownership, never
            // validity: fetch must still verify the recovered bytes before reuse.
            if !recorded.exists() && completed.is_file() && object.relative_path.starts_with("tmp/")
            {
                let o = self.index.objects.get_mut(&key).unwrap();
                o.relative_path = format!("{folder}/{key}");
                o.state = State::Partial;
                changed = true;
            } else if object.state == State::Downloading {
                self.index.objects.get_mut(&key).unwrap().state = State::Partial;
                changed = true;
            }
        }
        if changed {
            self.save_index()?;
        }
        // Exclusive ownership proves no cooperating workflow is active. A previous
        // running experiment is interrupted, even if some result files exist.
        for (path, _) in scan(&self.root.join("manifests"))? {
            if path.extension().is_some_and(|x| x == "json") {
                let mut m: Manifest = read_json(&path)?;
                m.validate()?;
                let mut interrupted = false;
                for e in &mut m.experiments {
                    if e.state == "running" {
                        e.state = "interrupted".into();
                        e.error = Some(
                            "Previous process exited before committing verified results".into(),
                        );
                        interrupted = true;
                    }
                }
                for attempt in &mut m.preparation_attempts {
                    if attempt.state == "running" {
                        attempt.state = "interrupted".into();
                        attempt.error =
                            Some("Previous process exited before committing prepared PIN".into());
                        interrupted = true;
                    }
                }
                if interrupted {
                    self.save_manifest(&m)?;
                }
            }
        }
        Ok(())
    }
    fn validate_index(&self) -> Result<()> {
        if self.index.schema_version != 1 {
            return Err("unsupported cache index schema".into());
        }
        let mut paths = BTreeSet::new();
        for (key, o) in &self.index.objects {
            if key != &o.key || !valid_key(key) || !paths.insert(&o.relative_path) {
                return Err("invalid/duplicate cache object identity".into());
            }
            let expected = [
                format!("objects/{key}"),
                format!("prepared/{key}"),
                format!("tmp/{key}.part"),
            ];
            if !expected.contains(&o.relative_path) {
                return Err("invalid cache index path; refusing filesystem access".into());
            }
            self.path(&o.relative_path)?;
        }
        Ok(())
    }
    pub fn path(&self, relative: &str) -> Result<PathBuf> {
        let p = Path::new(relative);
        if p.components().any(|c| !matches!(c, Component::Normal(_)))
            || p.as_os_str().is_empty()
            || relative.contains('\\')
        {
            return Err("unsafe managed-cache relative path".into());
        }
        let p = self.root.join(p);
        reject_symlinks(&p)?;
        Ok(p)
    }
    pub fn save_index(&self) -> Result<()> {
        self.require_write()?;
        self.validate_index()?;
        atomic_json(&self.path("index.json")?, &self.index)
    }
    pub fn require_write(&self) -> Result<()> {
        if self.writable {
            Ok(())
        } else {
            Err("read-only/dry-run cache cannot be modified".into())
        }
    }
    pub fn load_manifest(&self, a: &Pxd) -> Result<Manifest> {
        let m: Manifest = read_json(&self.path(&format!("manifests/{a}.json"))?)?;
        m.validate()?;
        Ok(m)
    }
    pub fn save_manifest(&self, m: &Manifest) -> Result<()> {
        self.require_write()?;
        m.validate()?;
        atomic_json(&self.path(&format!("manifests/{}.json", m.accession))?, m)
    }
    pub fn pin(&mut self, a: &Pxd, pin: bool) -> Result<()> {
        self.require_write()?;
        if pin {
            self.index.pinned.insert(a.clone());
        } else {
            self.index.pinned.remove(a);
        }
        self.save_index()
    }
    pub fn pinned(&self, o: &Object) -> bool {
        o.projects.iter().any(|p| self.index.pinned.contains(p))
    }
    pub fn evictable(&self, o: &Object, purge: bool) -> bool {
        !self.pinned(o)
            && o.reproducible
            && (purge
                || match o.retention {
                    Retention::Keep => false,
                    Retention::UntilResultVerified => o.result_verified,
                    Retention::Evict | Retention::KeepIfPinned => true,
                })
    }
    pub fn status(&self) -> Result<Status> {
        let sources = scan(&self.root.join("objects"))?;
        let prepared = scan(&self.root.join("prepared"))?;
        let temporary = scan(&self.root.join("tmp"))?;
        let all: Vec<_> = sources.iter().chain(&prepared).chain(&temporary).collect();
        let mut pinned = 0;
        let mut evictable = 0;
        let mut protected = 0;
        let mut tracked = BTreeSet::new();
        let mut largest = Vec::new();
        for (key, o) in &self.index.objects {
            let path = self.path(&o.relative_path)?;
            if path.is_file() {
                let len = fs::metadata(&path)?.len();
                tracked.insert(path);
                if self.pinned(o) {
                    pinned += len;
                } else if self.evictable(o, false) {
                    evictable += len;
                } else {
                    protected += len;
                }
                largest.push((key.clone(), len));
            }
        }
        largest.sort_by_key(|x| std::cmp::Reverse(x.1));
        largest.truncate(10);
        let sum = |items: &Vec<(PathBuf, u64)>| total(items.iter().map(|x| x.1));
        let source_bytes = sum(&sources)?;
        let prepared_bytes = sum(&prepared)?;
        let temporary_partial_bytes = sum(&temporary)?;
        let metadata = scan(&self.root.join("manifests"))?;
        Ok(Status {
            root: self.root.clone(),
            limit_bytes: self.limit,
            large_data_bytes: total([source_bytes, prepared_bytes, temporary_partial_bytes])?,
            source_bytes,
            prepared_bytes,
            temporary_partial_bytes,
            pinned_bytes: pinned,
            evictable_bytes: evictable,
            protected_by_retention_bytes: protected,
            untracked_bytes: total(all.iter().filter(|x| !tracked.contains(&x.0)).map(|x| x.1))?,
            metadata_bytes: sum(&metadata)?,
            results_bytes: sum(&scan(&self.root.join("results"))?)?,
            project_count: metadata
                .iter()
                .filter(|(p, _)| p.extension().is_some_and(|e| e == "json"))
                .count(),
            free_filesystem_bytes: available_space(&self.root)?,
            largest_objects: largest,
        })
    }
    pub fn eviction_plan(
        &self,
        additional: u64,
        free_needed: u64,
        exclude: &BTreeSet<String>,
    ) -> Result<Vec<String>> {
        let s = self.status()?;
        let capacity = total([s.large_data_bytes, additional])?.saturating_sub(self.limit);
        let disk = free_needed.saturating_sub(s.free_filesystem_bytes);
        let needed = capacity.max(disk);
        if needed == 0 {
            return Ok(vec![]);
        }
        let mut candidates: Vec<_> = self
            .index
            .objects
            .values()
            .filter(|o| !exclude.contains(&o.key) && self.evictable(o, false))
            .collect();
        candidates.sort_by_key(|o| (o.last_used_unix_seconds, &o.key));
        let mut out = Vec::new();
        let mut freed = 0;
        for o in candidates {
            let path = self.path(&o.relative_path)?;
            if path.is_file() {
                freed = total([freed, fs::metadata(path)?.len()])?;
                out.push(o.key.clone());
            }
            if freed >= needed {
                return Ok(out);
            }
        }
        Err(format!("storage preflight failed: need to free {needed} bytes; only {freed} safely evictable (limit {}, pinned {}, retained {}, untracked {})",self.limit,s.pinned_bytes,s.protected_by_retention_bytes,s.untracked_bytes).into())
    }
    pub fn evict(&mut self, keys: &[String], purge: bool) -> Result<()> {
        self.require_write()?;
        for key in keys {
            let o = self.index.objects.get(key).ok_or("unknown cache object")?;
            if !self.evictable(o, purge) {
                return Err(format!("{key} is pinned or protected by retention").into());
            }
            let path = self.path(&o.relative_path)?;
            // Record eviction before unlink. A crash can leave disposable bytes, never
            // a manifest falsely declaring a deleted artifact locally available.
            for (p, _) in scan(&self.root.join("manifests"))? {
                if p.extension().is_some_and(|x| x == "json") {
                    let mut m: Manifest = read_json(&p)?;
                    m.validate()?;
                    let mut changed = false;
                    for f in m.local_files.values_mut().filter(|f| &f.object_key == key) {
                        f.state = State::Evicted;
                        changed = true;
                    }
                    if changed {
                        self.save_manifest(&m)?;
                    }
                }
            }
            self.index.objects.get_mut(key).unwrap().state = State::Evicted;
            self.save_index()?;
            if path.exists() {
                fs::remove_file(path)?;
            }
        }
        Ok(())
    }
    /// Ordinary pruning drops disposable artifacts too: the ceiling is never a retention target.
    pub fn prune(&mut self, purge: bool, dry_run: bool, only_temporary: bool) -> Result<Cleanup> {
        let before = self.status()?;
        let mut keys = Vec::new();
        let mut freed = 0;
        for (key, o) in &self.index.objects {
            if self.evictable(o, purge) && (!only_temporary || o.relative_path.starts_with("tmp/"))
            {
                let path = self.path(&o.relative_path)?;
                if path.is_file() {
                    freed = total([freed, fs::metadata(path)?.len()])?;
                    keys.push(key.clone());
                }
            }
        }
        if !dry_run {
            self.evict(&keys, purge)?;
        }
        Ok(Cleanup {
            before_bytes: before.large_data_bytes,
            objects: keys,
            freed_bytes: freed,
            remaining_bytes: before.large_data_bytes - freed,
            pinned_remaining_bytes: before.pinned_bytes,
            provenance_results_retained: true,
            dry_run,
        })
    }
    pub fn record_download(
        &mut self,
        m: &mut Manifest,
        file: &RemoteFile,
        o: Object,
    ) -> Result<()> {
        m.local_files.insert(file.id.clone(), o.local());
        self.index.objects.insert(o.key.clone(), o);
        self.save_index()?;
        self.save_manifest(m)
    }
}
pub fn valid_key(key: &str) -> bool {
    let Some((kind, hash)) = key.split_once('-') else {
        return false;
    };
    matches!(kind, "md5" | "sha1" | "sha256" | "remote" | "prepared")
        && matches!(hash.len(), 32 | 40 | 64)
        && hash.bytes().all(|b| b.is_ascii_hexdigit())
}
pub fn reject_symlinks(path: &Path) -> Result<()> {
    let mut p = PathBuf::new();
    for part in path.components() {
        if matches!(part, Component::ParentDir | Component::CurDir) {
            return Err("cache path must not contain . or .. components".into());
        }
        p.push(part);
        match fs::symlink_metadata(&p) {
            Ok(m) if m.file_type().is_symlink() => {
                return Err(format!(
                    "symlinks are not allowed in managed-cache paths: {}",
                    p.display()
                )
                .into())
            }
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
    }
    Ok(())
}
pub fn available_space(path: &Path) -> Result<u64> {
    let existing = path
        .ancestors()
        .find(|p| p.exists())
        .ok_or("no existing filesystem ancestor")?;
    Ok(fs2::available_space(existing)?)
}
fn scan(dir: &Path) -> Result<Vec<(PathBuf, u64)>> {
    reject_symlinks(dir)?;
    if !dir.exists() {
        return Ok(vec![]);
    }
    let mut files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let kind = entry.file_type()?;
        if kind.is_symlink() {
            return Err(format!(
                "refusing symlink in managed cache: {}",
                entry.path().display()
            )
            .into());
        }
        if kind.is_dir() {
            files.extend(scan(&entry.path())?);
        } else if kind.is_file() {
            files.push((entry.path(), entry.metadata()?.len()));
        } else {
            return Err("special files are not allowed in the PRIDE cache".into());
        }
    }
    Ok(files)
}
pub fn read_json<T: serde::de::DeserializeOwned>(p: &Path) -> Result<T> {
    reject_symlinks(p)?;
    if fs::metadata(p)?.len() > 64 * 1024 * 1024 {
        return Err("cache metadata exceeds 64 MiB safety limit".into());
    }
    Ok(serde_json::from_reader(File::open(p)?)?)
}
pub fn atomic_json<T: Serialize>(p: &Path, v: &T) -> Result<()> {
    reject_symlinks(p)?;
    let tmp = p.with_extension("json.new");
    reject_symlinks(&tmp)?;
    let mut f = OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(&tmp)?;
    serde_json::to_writer_pretty(&mut f, v)?;
    f.write_all(b"\n")?;
    f.sync_all()?;
    fs::rename(tmp, p)?;
    File::open(p.parent().ok_or("metadata parent missing")?)?.sync_all()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn releasing_cache_unlocks_an_inherited_file_description() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("cache");
        let cache = Cache::open(&root, 100, true).unwrap();
        let inherited = cache._lock.as_ref().unwrap().try_clone().unwrap();
        drop(cache);
        let reopened = Cache::open(&root, 100, true).unwrap();
        drop(reopened);
        drop(inherited);
    }
}
