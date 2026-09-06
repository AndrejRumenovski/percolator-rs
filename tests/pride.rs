use percolator_rs::pride::{
    self,
    cache::{Cache, Object, DEFAULT_LIMIT},
    client::{self, PrideClient},
    download::{self, Budgets, Downloader},
    workflow::{self, RunOptions},
    *,
};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs,
    io::{Read, Write},
    net::TcpListener,
    path::Path,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
    thread,
    time::Duration,
};
use tempfile::TempDir;

struct Reply {
    status: u16,
    body: Vec<u8>,
    length: Option<usize>,
    headers: Vec<(String, String)>,
}
impl Reply {
    fn ok(body: impl Into<Vec<u8>>) -> Self {
        Self {
            status: 200,
            body: body.into(),
            length: None,
            headers: vec![],
        }
    }
}
struct Server {
    url: String,
    requests: Arc<Mutex<Vec<String>>>,
    stop: Arc<AtomicBool>,
    thread: Option<thread::JoinHandle<()>>,
}
impl Server {
    fn new(handler: impl Fn(&str) -> Reply + Send + 'static) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let url = format!("http://{}", listener.local_addr().unwrap());
        let stop = Arc::new(AtomicBool::new(false));
        let done = stop.clone();
        let requests = Arc::new(Mutex::new(Vec::new()));
        let seen = requests.clone();
        let thread = thread::spawn(move || {
            while !done.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        stream
                            .set_read_timeout(Some(Duration::from_secs(2)))
                            .unwrap();
                        let mut request = Vec::new();
                        let mut buf = [0u8; 2048];
                        while let Ok(n) = stream.read(&mut buf) {
                            if n == 0 {
                                break;
                            }
                            request.extend_from_slice(&buf[..n]);
                            if request.windows(4).any(|s| s == b"\r\n\r\n") {
                                break;
                            }
                        }
                        let request = String::from_utf8_lossy(&request).into_owned();
                        seen.lock().unwrap().push(request.clone());
                        let reply = handler(&request);
                        let mut header = format!(
                            "HTTP/1.1 {} test\r\nConnection: close\r\nContent-Length: {}\r\n",
                            reply.status,
                            reply.length.unwrap_or(reply.body.len())
                        );
                        for (k, v) in reply.headers {
                            header.push_str(&format!("{k}: {v}\r\n"));
                        }
                        header.push_str("\r\n");
                        let _ = stream.write_all(header.as_bytes());
                        for chunk in reply.body.chunks(997) {
                            if stream.write_all(chunk).is_err() {
                                break;
                            }
                        }
                        let _ = stream.flush();
                    }
                    Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(2))
                    }
                    Err(_) => break,
                }
            }
        });
        Self {
            url,
            requests,
            stop,
            thread: Some(thread),
        }
    }
}
impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        self.thread.take().unwrap().join().unwrap();
    }
}
fn pxd() -> Pxd {
    "PXD000001".parse().unwrap()
}
fn manifest() -> Manifest {
    Manifest::new(
        client::parse_project(serde_json::json!({"accession":"PXD000001"})).unwrap(),
        vec![],
        API_BASE.into(),
    )
}
fn budgets() -> Budgets {
    Budgets {
        max_download: 10_000_000,
        max_working_space: None,
        safety: 0,
    }
}
fn cache(temp: &TempDir, limit: u64) -> Cache {
    Cache::open(&temp.path().join("relocated/cache"), limit, true).unwrap()
}
fn file(url: &str, body: &[u8]) -> RemoteFile {
    RemoteFile {
        id: "remote-id".into(),
        filename: "sample.pin".into(),
        category: Some("RESULT".into()),
        format: None,
        size_bytes: Some(body.len() as u64),
        checksum_table_size: None,
        references: vec![url.into()],
        checksums: vec![Checksum::new(
            Algorithm::Sha256,
            &format!("{:x}", Sha256::digest(body)),
            "mock repository",
        )
        .unwrap()],
        untyped_checksum: None,
        analysis_accessions: None,
        run_metadata: None,
        inventory_source: "fixture".into(),
    }
}
fn downloader() -> Downloader {
    Downloader::new(&budgets(), Arc::new(AtomicBool::new(false))).unwrap()
}
fn put(c: &mut Cache, m: &mut Manifest, body: &[u8], age: u64, retention: Retention) -> String {
    let f = file("https://ftp.pride.ebi.ac.uk/test", body);
    let key = f.object_key();
    let rel = format!("objects/{key}");
    fs::write(c.path(&rel).unwrap(), body).unwrap();
    let o = Object {
        key: key.clone(),
        relative_path: rel,
        bytes: body.len() as u64,
        state: State::Verified,
        local_sha256: Some(format!("{:x}", Sha256::digest(body))),
        verification: vec![],
        projects: BTreeSet::from([m.accession.clone()]),
        last_used_unix_seconds: age,
        retention,
        result_verified: false,
        reproducible: true,
        etag: None,
        pin_validated: false,
    };
    c.record_download(m, &f, o).unwrap();
    key
}

#[test]
fn accession_validation_and_serde_cannot_bypass_it() {
    for s in ["PXD000001", "PXD012345", "PXD999999"] {
        assert_eq!(s.parse::<Pxd>().unwrap().to_string(), s);
    }
    for s in [
        "",
        "PXD1",
        "pxd000001",
        "PXD0000012",
        "PXD１２３４５６",
        "PXD00000/",
        "/PXD000001",
        "PXD000001 ",
        "PXD../001",
    ] {
        assert!(s.parse::<Pxd>().is_err(), "{s}");
        assert!(serde_json::from_value::<Pxd>(serde_json::json!(s)).is_err());
    }
}
#[test]
fn missing_metadata_remains_null_and_bad_shapes_fail() {
    let m = manifest();
    assert!(m.project.title.is_none());
    assert!(m.project.organisms.is_none());
    assert_eq!(
        serde_json::to_value(&m.project).unwrap()["organisms"],
        serde_json::Value::Null
    );
    assert!(client::parse_project(
        serde_json::json!({"accession":"PXD000001","organisms":"human"})
    )
    .is_err());
    assert!(client::parse_files(
        serde_json::json!([{"accession":"id","fileName":"a","fileSizeBytes":-1}])
    )
    .is_err());
}
#[test]
fn frozen_official_metadata_checksum_table_and_inventory() {
    let p = client::parse_project(
        serde_json::from_str(include_str!("fixtures/pride/project.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(
        p.organisms.as_ref().unwrap()[0].name.as_deref(),
        Some("Erwinia carotovora")
    );
    let files = client::parse_files(
        serde_json::from_str(include_str!("fixtures/pride/files.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(files.len(), 8);
    let mut m = Manifest::new(p, files, API_BASE.into());
    client::merge_checksums(
        &mut m,
        include_str!("fixtures/pride/checksums.tsv"),
        "official MD5 endpoint",
        Some("ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001"),
    )
    .unwrap();
    assert_eq!(m.inventory.len(), 13);
    assert_eq!(m.indexed_file_count, 8);
    let readme = m
        .inventory
        .iter()
        .find(|f| f.filename == "README.txt")
        .unwrap();
    assert!(readme.references[0].ends_with("/PXD000001/README.txt"));
    let fasta = m
        .inventory
        .iter()
        .find(|f| f.filename.ends_with(".fasta"))
        .unwrap();
    assert_eq!(fasta.checksums[0].value.len(), 32);
    assert_eq!(fasta.checksums[0].reported_value.len(), 32);
    let conflict = m
        .inventory
        .iter()
        .find(|f| f.filename == "PRIDE_Exp_Complete_Ac_22134.xml.gz")
        .unwrap();
    assert!(conflict
        .size()
        .unwrap_err()
        .to_string()
        .contains("contradictory"));
    assert!(!m
        .inventory
        .iter()
        .any(|f| f.compatibility() == Compatibility::DirectlyCompatible));
}
#[test]
fn api_pagination_follows_count_even_if_server_caps_page_size() {
    let server = Server::new(|r| {
        let path = r.split_whitespace().nth(1).unwrap();
        if path == "/projects/PXD000001/files/count" {
            Reply::ok("3")
        } else {
            let page = path.split("page=").nth(1).unwrap();
            Reply::ok(format!(
                "[{{\"accession\":\"{page}\",\"fileName\":\"f{page}.raw\"}}]"
            ))
        }
    });
    let files = PrideClient::with_base(&server.url)
        .unwrap()
        .files(&pxd())
        .unwrap();
    assert_eq!(files.len(), 3);
    let requests = server.requests.lock().unwrap();
    assert_eq!(requests.len(), 5);
    assert!(requests[3].contains("page=2"));
}
#[test]
fn api_repeated_page_and_early_empty_are_errors() {
    for body in ["[]", "[{\"accession\":\"same\",\"fileName\":\"x.pin\"}]"] {
        let server = Server::new(move |r| {
            Reply::ok(if r.contains("/files/count ") {
                "3"
            } else {
                body
            })
        });
        assert!(PrideClient::with_base(&server.url)
            .unwrap()
            .files(&pxd())
            .is_err());
    }
}
#[test]
fn api_retry_status_search_pagination_and_error_diagnostics() {
    let calls = Arc::new(Mutex::new(0));
    let counter = calls.clone();
    let server = Server::new(move |r| {
        if r.contains("/status/") {
            Reply::ok("PRIVATE")
        } else {
            let mut n = counter.lock().unwrap();
            *n += 1;
            if *n == 1 {
                Reply {
                    status: 503,
                    ..Reply::ok("busy")
                }
            } else {
                Reply::ok("[{\"accession\":\"PXD000002\",\"organisms\":[\"Homo sapiens\"]}]")
            }
        }
    });
    let client = PrideClient::with_base(&server.url).unwrap();
    assert!(client
        .project(&pxd())
        .unwrap_err()
        .to_string()
        .contains("PRIVATE"));
    let r = client
        .search("heart & brain", Some("organisms==Homo sapiens"), 2, 1)
        .unwrap();
    assert_eq!(r.len(), 1);
    assert_eq!(*calls.lock().unwrap(), 2);
    assert!(server
        .requests
        .lock()
        .unwrap()
        .iter()
        .any(|s| s.contains("page=2") && s.contains("heart+%26+brain")));
}
#[test]
fn filters_and_compatibility_are_conservative() {
    let mut f = file("https://example.org/x", b"x");
    for (name, kind) in [
        ("x.pin", Compatibility::PotentiallyConvertible),
        ("x.pin.gz", Compatibility::PotentiallyConvertible),
        ("x.mzML", Compatibility::RawRequiresSearch),
        ("x.mzid", Compatibility::PotentiallyConvertible),
        ("x.mzTab", Compatibility::PotentiallyConvertible),
        ("x.raw", Compatibility::RawRequiresSearch),
        ("x.fasta", Compatibility::UnrelatedUnknown),
    ] {
        f.filename = name.into();
        assert_eq!(f.compatibility(), kind, "{name}");
    }
    f.filename = "run.pep.xml".into();
    assert!(f.matches("search-engine-output"));
    assert!(f.matches("processed"));
    f.filename = "x.mzid.gz".into();
    assert!(f.matches("mzIdentML"));
}
#[test]
fn complete_streaming_download_and_authoritative_verification() {
    let body = vec![73u8; 250_000];
    let sent = body.clone();
    let server = Server::new(move |_| Reply::ok(sent.clone()));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 300_000);
    let mut m = manifest();
    let f = file(&server.url, &body);
    let key = downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap();
    let o = &c.index.objects[&key];
    assert_eq!(o.state, State::Verified);
    assert!(o.verification[0].verified);
    assert_eq!(fs::read(c.path(&o.relative_path).unwrap()).unwrap(), body);
    assert_eq!(c.status().unwrap().large_data_bytes, 250_000);
    assert_eq!(c.status().unwrap().temporary_partial_bytes, 0);
}
#[test]
fn checksum_mismatch_never_publishes_verified_object() {
    let server = Server::new(|_| Reply::ok("BAD"));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 100);
    let mut m = manifest();
    let f = file(&server.url, b"abc");
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .is_err());
    let o = &c.index.objects[&f.object_key()];
    assert_eq!(o.state, State::Corrupt);
    assert!(!o.verification[0].verified);
    assert!(o.relative_path.ends_with(".part"));
    assert_eq!(m.local_files[&f.id].state, State::Corrupt);
    assert!(o.local_sha256.is_some());
}
#[test]
fn interrupted_download_resumes_only_valid_content_range() {
    let server = Server::new(|r| {
        if r.to_ascii_lowercase().contains("range: bytes=3-") {
            Reply {
                status: 206,
                headers: vec![("Content-Range".into(), "bytes 3-5/6".into())],
                ..Reply::ok("def")
            }
        } else {
            Reply {
                length: Some(6),
                ..Reply::ok("abc")
            }
        }
    });
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let f = file(&server.url, b"abcdef");
    let key = downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap();
    assert_eq!(c.index.objects[&key].state, State::Verified);
    assert_eq!(server.requests.lock().unwrap().len(), 2);
}
#[test]
fn invalid_range_and_cancelled_download_remain_partial() {
    let server = Server::new(|r| {
        if r.to_ascii_lowercase().contains("range:") {
            Reply {
                status: 206,
                headers: vec![("Content-Range".into(), "bytes 0-2/6".into())],
                ..Reply::ok("def")
            }
        } else {
            Reply {
                length: Some(6),
                ..Reply::ok("abc")
            }
        }
    });
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let f = file(&server.url, b"abcdef");
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap_err()
        .to_string()
        .contains("Content-Range"));
    assert_eq!(c.index.objects[&f.object_key()].state, State::Partial);
    let mut dl = downloader();
    dl.cancelled.store(true, Ordering::Relaxed);
    assert!(dl.fetch(&mut c, &mut m, &f, &BTreeSet::new()).is_err());
    let cleanup = c.prune(false, false, true).unwrap();
    assert_eq!(cleanup.freed_bytes, 3);
    assert_eq!(cleanup.remaining_bytes, 0);
}
#[test]
fn ignored_range_restarts_and_missing_checksum_is_explicit() {
    let calls = Arc::new(Mutex::new(0));
    let n = calls.clone();
    let server = Server::new(move |_| {
        let mut n = n.lock().unwrap();
        *n += 1;
        if *n == 1 {
            Reply {
                length: Some(6),
                headers: vec![("ETag".into(), "\"version1\"".into())],
                ..Reply::ok("abc")
            }
        } else {
            Reply::ok("abcdef")
        }
    });
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let mut f = file(&server.url, b"abcdef");
    f.checksums.clear();
    let key = downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap();
    assert_eq!(c.index.objects[&key].state, State::DownloadedUnverified);
    assert!(m.local_files[&f.id].verification_unavailable);
    assert!(server.requests.lock().unwrap()[1]
        .to_ascii_lowercase()
        .contains("if-range: \"version1\""));
}
#[test]
fn content_dedup_cross_project_and_reverification() {
    let server = Server::new(|_| Reply::ok("abc"));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let f = file(&server.url, b"abc");
    let key = downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap();
    let mut other = manifest();
    other.accession = "PXD000002".parse().unwrap();
    other.project.accession = other.accession.clone();
    let mut alias = f.clone();
    alias.id = "different-record".into();
    alias.filename = "different.pin".into();
    assert_eq!(
        downloader()
            .fetch(&mut c, &mut other, &alias, &BTreeSet::new())
            .unwrap(),
        key
    );
    assert_eq!(server.requests.lock().unwrap().len(), 1);
    assert_eq!(c.status().unwrap().large_data_bytes, 3);
    assert_eq!(c.index.objects[&key].projects.len(), 2);
    fs::write(
        c.path(&c.index.objects[&key].relative_path).unwrap(),
        b"bad",
    )
    .unwrap();
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .is_err());
    assert_eq!(c.index.objects[&key].state, State::Corrupt);
}
#[test]
fn lru_eviction_pins_and_hard_ceiling() {
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let old = put(&mut c, &mut m, b"aaaa", 1, Retention::Evict);
    let new = put(&mut c, &mut m, b"bbbb", 2, Retention::Evict);
    let plan = c.eviction_plan(5, 0, &BTreeSet::new()).unwrap();
    assert_eq!(plan, vec![old.clone()]);
    c.evict(&plan, false).unwrap();
    assert!(c.index.objects.contains_key(&old));
    assert_eq!(c.index.objects[&old].state, State::Evicted);
    c.pin(&pxd(), true).unwrap();
    assert_eq!(c.status().unwrap().pinned_bytes, 4);
    assert!(c.eviction_plan(8, 0, &BTreeSet::new()).is_err());
    assert!(c.evict(&[new], true).is_err());
    c.pin(&pxd(), false).unwrap();
    assert_eq!(c.prune(false, false, false).unwrap().remaining_bytes, 0);
}
#[test]
fn download_working_space_free_space_and_default_limit() {
    assert_eq!(DEFAULT_LIMIT, 50_000_000_000);
    assert_eq!(pride::bytes("50GB").unwrap(), DEFAULT_LIMIT);
    assert!(pride::bytes("1.5GB").is_err());
    assert!(pride::bytes("18446744073709551615TB").is_err());
    let tmp = TempDir::new().unwrap();
    let c = cache(&tmp, 10);
    let m = manifest();
    let f = file("https://example.org/x", b"12345678");
    let mut b = budgets();
    b.max_download = 7;
    assert!(
        download::plan(&c, &m, std::slice::from_ref(&f), &b, false, 1, 0)
            .unwrap_err()
            .to_string()
            .contains("download budget")
    );
    b.max_download = 8;
    b.max_working_space = Some(9);
    assert!(
        download::plan(&c, &m, std::slice::from_ref(&f), &b, false, 1, 2)
            .unwrap_err()
            .to_string()
            .contains("working-space")
    );
    b.max_working_space = None;
    assert!(download::plan(
        &c,
        &m,
        &[file("https://example.org/x", b"12345678901")],
        &b,
        false,
        1,
        0
    )
    .is_err());
    assert!(c
        .eviction_plan(0, u64::MAX, &BTreeSet::new())
        .unwrap_err()
        .to_string()
        .contains("storage preflight"));
    let before = fs::read_dir(&c.root).unwrap().count();
    let p = download::plan(&c, &m, &[f], &b, false, 1, 0).unwrap();
    assert_eq!(p.expected_final_large_data_bytes, 8);
    assert_eq!(fs::read_dir(&c.root).unwrap().count(), before);
}
#[test]
fn cleanup_retention_purge_and_manifest_lineage_survive() {
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 100);
    let mut m = manifest();
    let keep = put(&mut c, &mut m, b"keep", 1, Retention::Keep);
    let until = put(&mut c, &mut m, b"until", 2, Retention::UntilResultVerified);
    put(&mut c, &mut m, b"discard", 3, Retention::Evict);
    m.lineage.push(Lineage {
        id: "pin".into(),
        inputs: vec!["raw".into()],
        output_sha256: Some("hash".into()),
        kind: "pin".into(),
        tool: "external converter".into(),
        tool_version: Some("1.0".into()),
        parameters: vec!["--explicit-decoys".into()],
        protein_database: Some("db.fasta".into()),
        database_sha256: Some("database hash".into()),
        decoy_generation: Some("recipe".into()),
    });
    c.save_manifest(&m).unwrap();
    fs::write(c.path("results/final.tsv").unwrap(), b"final").unwrap();
    let preview = c.prune(false, true, false).unwrap();
    assert_eq!(preview.freed_bytes, 7);
    assert_eq!(c.status().unwrap().large_data_bytes, 16);
    let report = c.prune(false, false, false).unwrap();
    assert_eq!(report.remaining_bytes, 9);
    assert_eq!(
        c.load_manifest(&pxd()).unwrap().lineage[0].inputs,
        vec!["raw"]
    );
    c.index.objects.get_mut(&until).unwrap().result_verified = true;
    c.save_index().unwrap();
    assert_eq!(c.prune(false, false, false).unwrap().remaining_bytes, 4);
    assert_eq!(c.prune(true, false, false).unwrap().remaining_bytes, 0);
    assert_eq!(c.index.objects[&keep].state, State::Evicted);
    assert_eq!(
        fs::read(c.path("results/final.tsv").unwrap()).unwrap(),
        b"final"
    );
    assert_eq!(c.load_manifest(&pxd()).unwrap().lineage.len(), 1);
    let json = serde_json::to_vec(&c.load_manifest(&pxd()).unwrap()).unwrap();
    let round: Manifest = serde_json::from_slice(&json).unwrap();
    round.validate().unwrap();
}
#[test]
fn ephemeral_planning_allows_project_larger_than_cache_but_pins_must_fit() {
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let m = manifest();
    let files = vec![
        file("https://example.org/a", b"aaaaaaa"),
        file("https://example.org/b", b"bbbbbbb"),
    ];
    assert!(download::plan(&c, &m, &files, &budgets(), false, 1, 0).is_err());
    let plan = download::plan(&c, &m, &files, &budgets(), true, 1, 0).unwrap();
    assert_eq!(plan.download_bytes, 14);
    assert_eq!(plan.expected_final_large_data_bytes, 0);
    c.pin(&pxd(), true).unwrap();
    assert!(download::plan(&c, &m, &files, &budgets(), true, 1, 0).is_err());
}
#[test]
fn dry_cache_does_not_create_any_files_and_relocation_is_explicit() {
    let tmp = TempDir::new().unwrap();
    let root = tmp.path().join("external/cache");
    let c = Cache::open(&root, DEFAULT_LIMIT, false).unwrap();
    assert_eq!(c.status().unwrap().large_data_bytes, 0);
    assert!(!root.exists());
    assert!(c.save_index().is_err());
    let owned = cache(&tmp, 10);
    assert!(owned.root.starts_with(tmp.path()));
    assert!(Cache::open(&owned.root, 10, true).is_err());
}
#[test]
fn unsafe_remote_names_cannot_escape_and_cache_symlinks_are_rejected() {
    for name in ["../a", "/a", "x/../a", "C:\\evil", "a\\b", "a\n"] {
        assert!(!client::safe_remote_relative(name));
    }
    let server = Server::new(|_| Reply::ok("abc"));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let mut f = file(&server.url, b"abc");
    f.filename = "../../outside.pin".into();
    let key = downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .unwrap();
    assert!(c.index.objects[&key]
        .relative_path
        .starts_with("objects/sha256-"));
    assert!(!tmp.path().join("outside.pin").exists());
    assert!(c.path("../outside").is_err());
    assert!(c.path("/tmp/outside").is_err());
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(tmp.path(), c.root.join("tmp/escape")).unwrap();
        assert!(c.prune(true, false, false).is_err());
        assert!(c.path("tmp/escape/anything").is_err());
    }
}
#[test]
fn corrupt_index_and_unowned_directory_fail_closed() {
    let tmp = TempDir::new().unwrap();
    let c = cache(&tmp, 10);
    let root = c.root.clone();
    drop(c);
    fs::write(root.join("index.json"), b"{broken").unwrap();
    let error = Cache::open(&root, 10, true).unwrap_err_text();
    assert!(error.contains("corrupt"), "unexpected error: {error}");
    let unowned = tmp.path().join("unowned");
    fs::create_dir(&unowned).unwrap();
    fs::write(unowned.join("keep"), b"x").unwrap();
    assert!(Cache::open(&unowned, 10, true).is_err());
    assert!(unowned.join("keep").exists());
}
trait ErrorText {
    fn unwrap_err_text(self) -> String;
}
impl<T> ErrorText for Result<T> {
    fn unwrap_err_text(self) -> String {
        match self {
            Ok(_) => panic!("expected failure"),
            Err(e) => e.to_string(),
        }
    }
}
#[test]
fn rename_crash_is_recovered_as_owned_partial_and_prunable() {
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 10);
    let mut m = manifest();
    let key = put(&mut c, &mut m, b"abc", 0, Retention::Evict);
    c.index.objects.get_mut(&key).unwrap().relative_path = format!("tmp/{key}.part");
    c.index.objects.get_mut(&key).unwrap().state = State::Downloading;
    c.save_index().unwrap();
    let root = c.root.clone();
    drop(c);
    let mut c = Cache::open(&root, 10, true).unwrap();
    assert_eq!(c.index.objects[&key].state, State::Partial);
    assert_eq!(c.prune(false, false, false).unwrap().remaining_bytes, 0);
}
#[test]
fn scientific_argument_scope_and_raw_run_refusals() {
    for args in [
        vec!["--join"],
        vec!["--results-psms", "/tmp/x"],
        vec!["other.pin"],
        vec!["--maxiter", "bad"],
        vec!["--profile", "bad"],
    ] {
        assert!(workflow::validate_analysis_args(
            &args.into_iter().map(str::to_owned).collect::<Vec<_>>()
        )
        .is_err());
    }
    workflow::validate_analysis_args(&[
        "--profile".into(),
        "fast".into(),
        "--seed".into(),
        "5".into(),
    ])
    .unwrap();
    let mut f = file("https://example.org/a", b"pin");
    let mut options = RunOptions::default();
    assert!(workflow::validate_run(&[f.clone(), f.clone()], &options).is_err());
    options.independent_runs = true;
    f.filename = "spectra.mzML".into();
    assert!(workflow::validate_run(&[f], &options).is_err());
}

#[test]
fn real_analysis_ephemeral_outputs_equal_existing_cli_and_sources_reach_zero() {
    let body = include_bytes!("fixtures/sample.pin").to_vec();
    let sent = body.clone();
    let server = Server::new(move |_| Reply::ok(sent.clone()));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, body.len() as u64 + 64 * 1024 * 1024);
    let mut m = manifest();
    let f = file(&server.url, &body);
    m.inventory = vec![f.clone()];
    c.save_manifest(&m).unwrap();
    let opts = RunOptions {
        ephemeral: true,
        analysis_args: vec!["--maxiter".into(), "1".into()],
        ..RunOptions::default()
    };
    workflow::run(
        &mut c,
        &mut m,
        std::slice::from_ref(&f),
        &mut downloader(),
        &budgets(),
        &opts,
        Path::new(env!("CARGO_BIN_EXE_percolator-rs")),
    )
    .unwrap();
    assert_eq!(c.status().unwrap().large_data_bytes, 0);
    let m = c.load_manifest(&pxd()).unwrap();
    let e = &m.experiments[0];
    assert_eq!(e.state, "verified");
    assert_eq!(e.result_hashes.len(), 4);
    assert_eq!(e.lineage.len(), 5);
    assert!(e.executable_sha256.len() == 64);
    assert_eq!(m.local_files[&f.id].state, State::Evicted);
    let baseline = tmp.path().join("baseline.tsv");
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_percolator-rs"))
        .args(["--maxiter", "1", "--results-psms"])
        .arg(&baseline)
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/sample.pin"
        ))
        .output()
        .unwrap();
    assert!(output.status.success());
    let result = e
        .result_hashes
        .keys()
        .find(|k| k.ends_with("/target.psms.tsv"))
        .unwrap();
    assert_eq!(
        fs::read(c.path(result).unwrap()).unwrap(),
        fs::read(baseline).unwrap()
    );
    c.prune(true, false, false).unwrap();
    assert_eq!(
        c.load_manifest(&pxd()).unwrap().experiments[0]
            .result_hashes
            .len(),
        4
    );
}
#[test]
fn analysis_failure_does_not_mark_complete_or_remove_verified_source() {
    let server = Server::new(|_| Reply::ok("not a pin\n"));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 64 * 1024 * 1024 + 100);
    let mut m = manifest();
    let f = file(&server.url, b"not a pin\n");
    assert!(workflow::run(
        &mut c,
        &mut m,
        &[f],
        &mut downloader(),
        &budgets(),
        &RunOptions {
            ephemeral: true,
            ..RunOptions::default()
        },
        Path::new(env!("CARGO_BIN_EXE_percolator-rs"))
    )
    .is_err());
    let m = c.load_manifest(&pxd()).unwrap();
    assert_eq!(m.experiments[0].state, "failed");
    assert!(m.experiments[0].result_hashes.is_empty());
    assert_eq!(c.status().unwrap().source_bytes, 10);
}
#[test]
#[ignore = "optional official PRIDE network smoke test; normal CI stays offline"]
fn live_official_pride_metadata() {
    let m = PrideClient::new().unwrap().manifest(&pxd()).unwrap();
    assert_eq!(m.project.status.as_deref(), Some("PUBLIC"));
    assert!(m.indexed_file_count > 0);
    assert!(m.inventory.len() >= m.indexed_file_count);
}

fn recipe(source: &str, id: &str) -> pride::prepare::Recipe {
    pride::prepare::Recipe {
        steps: vec![Lineage {
            id: id.into(),
            inputs: vec![source.into()],
            output_sha256: None,
            kind: "pin".into(),
            tool: "fixture exporter".into(),
            tool_version: Some("1.0".into()),
            parameters: vec![
                "export complete target/decoy candidates with numeric features".into(),
            ],
            protein_database: None,
            database_sha256: None,
            decoy_generation: None,
        }],
    }
}
#[test]
fn prepared_import_lineage_analysis_and_purge_are_real_operations() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("external.pin");
    fs::write(&input, include_bytes!("fixtures/sample.pin")).unwrap();
    let mut c = cache(&tmp, 100_000_000);
    let mut m = manifest();
    m.inventory
        .push(file("https://example.org/search-result", b"source"));
    c.save_manifest(&m).unwrap();
    let id = pride::prepare::import(
        &mut c,
        &mut m,
        &input,
        recipe("remote-id", "exported-pin"),
        &budgets(),
        Retention::KeepIfPinned,
    )
    .unwrap();
    assert_eq!(
        c.status().unwrap().prepared_bytes,
        fs::metadata(&input).unwrap().len()
    );
    assert_eq!(m.lineage[0].inputs, vec!["remote-id"]);
    let opts = RunOptions {
        ephemeral: true,
        analysis_args: vec!["--maxiter".into(), "1".into()],
        ..RunOptions::default()
    };
    workflow::run_prepared(
        &mut c,
        &mut m,
        &id,
        &Arc::new(AtomicBool::new(false)),
        &budgets(),
        &opts,
        Path::new(env!("CARGO_BIN_EXE_percolator-rs")),
    )
    .unwrap();
    assert_eq!(c.status().unwrap().large_data_bytes, 0);
    assert!(input.exists());
    assert_eq!(m.experiments[0].lineage[0].inputs, vec!["exported-pin"]);
    assert!(pride::prepare::run_plan(&c, &m, &id, &budgets(), 4096).is_err());
    c.prune(true, false, false).unwrap();
    let m = c.load_manifest(&pxd()).unwrap();
    assert_eq!(m.prepared_pins.len(), 1);
    assert_eq!(m.lineage.len(), 1);
    assert_eq!(m.experiments[0].state, "verified");
}
#[test]
fn preparation_requires_recipe_and_search_database_identity() {
    let mut m = manifest();
    m.inventory
        .push(file("https://example.org/source", b"source"));
    let mut r = recipe("wrong-id", "pin");
    assert!(pride::prepare::validate_recipe(&m, &r).is_err());
    r = recipe("remote-id", "pin");
    r.steps[0].tool_version = None;
    assert!(pride::prepare::validate_recipe(&m, &r).is_err());
    r = recipe("remote-id", "search");
    r.steps[0].kind = "database_search".into();
    r.steps.push(recipe("search", "pin").steps.remove(0));
    assert!(pride::prepare::validate_recipe(&m, &r).is_err());
    r.steps[0].protein_database = Some("database.fasta".into());
    r.steps[0].database_sha256 = Some("a".repeat(64));
    r.steps[0].decoy_generation = Some("reverse".into());
    pride::prepare::validate_recipe(&m, &r).unwrap();
}
#[test]
fn output_ceiling_failure_leaves_prunable_temporary_data_and_valid_source() {
    let body = include_bytes!("fixtures/sample.pin").to_vec();
    let sent = body.clone();
    let server = Server::new(move |_| Reply::ok(sent.clone()));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, body.len() as u64 + 4096);
    let mut m = manifest();
    let f = file(&server.url, &body);
    let opts = RunOptions {
        ephemeral: true,
        result_bytes_per_input: 4096,
        analysis_args: vec!["--maxiter".into(), "1".into()],
        ..RunOptions::default()
    };
    assert!(workflow::run(
        &mut c,
        &mut m,
        &[f],
        &mut downloader(),
        &budgets(),
        &opts,
        Path::new(env!("CARGO_BIN_EXE_percolator-rs"))
    )
    .is_err());
    assert_eq!(m.experiments[0].state, "failed");
    assert!(m.experiments[0].result_hashes.is_empty());
    let s = c.status().unwrap();
    assert!(s.large_data_bytes <= c.limit);
    assert_eq!(s.results_bytes, 0);
    assert!(s.temporary_partial_bytes <= 4096);
    assert_eq!(c.prune(false, false, false).unwrap().remaining_bytes, 0);
}
#[test]
fn independent_batches_process_more_source_bytes_than_ceiling() {
    let body = include_bytes!("fixtures/sample.pin").to_vec();
    let mut other = body.clone();
    other.extend_from_slice(b"\n");
    let a = body.clone();
    let b = other.clone();
    let server = Server::new(move |r| {
        if r.starts_with("GET /b ") {
            Reply::ok(b.clone())
        } else {
            Reply::ok(a.clone())
        }
    });
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, other.len() as u64 + 2_000_000);
    let mut m = manifest();
    let f = file(&format!("{}/a", server.url), &body);
    let mut g = file(&format!("{}/b", server.url), &other);
    g.id = "other-run".into();
    m.inventory = vec![f.clone(), g.clone()];
    let opts = RunOptions {
        ephemeral: true,
        independent_runs: true,
        batch_size: 1,
        result_bytes_per_input: 2_000_000,
        analysis_args: vec!["--maxiter".into(), "1".into()],
        ..RunOptions::default()
    };
    assert!(body.len() as u64 + other.len() as u64 > c.limit);
    workflow::run(
        &mut c,
        &mut m,
        &[f, g],
        &mut downloader(),
        &budgets(),
        &opts,
        Path::new(env!("CARGO_BIN_EXE_percolator-rs")),
    )
    .unwrap();
    assert_eq!(m.experiments.len(), 2);
    assert_eq!(c.status().unwrap().large_data_bytes, 0);
}
#[test]
fn oversized_response_and_download_ceiling_stop_before_writing_excess() {
    let server = Server::new(|_| Reply::ok("abcdef"));
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 3);
    let mut m = manifest();
    let f = file(&server.url, b"abc");
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .is_err());
    assert!(c.status().unwrap().large_data_bytes <= 3);
    let f = file(&server.url, b"abcdef");
    let count = server.requests.lock().unwrap().len();
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .is_err());
    assert_eq!(server.requests.lock().unwrap().len(), count);
}
#[test]
fn untyped_checksum_and_ambiguous_table_names_fail_closed() {
    assert!(Checksum::new(Algorithm::Md5, "abc", "PRIDE").is_err());
    assert!(Checksum::new(Algorithm::Sha1, &"z".repeat(40), "PRIDE").is_err());
    let mut m = manifest();
    let f = file("https://example.org/a", b"a");
    let mut g = f.clone();
    g.id = "second".into();
    m.inventory = vec![f, g];
    assert!(client::merge_checksums(
        &mut m,
        "File-Name\tFile-MD5Checksum\tFile-Size\nsample.pin\t0cc175b9c0f1b6a831c399e269772661\t1\n",
        "PRIDE",
        None
    )
    .is_err());
}
#[test]
fn cli_dry_run_prune_all_and_purge_confirmation() {
    let tmp = TempDir::new().unwrap();
    let root = tmp.path().join("cli-cache");
    let cmd = |args: &[&str]| {
        std::process::Command::new(env!("CARGO_BIN_EXE_percolator-rs"))
            .args(["pride", "--cache-dir"])
            .arg(&root)
            .args(args)
            .output()
            .unwrap()
    };
    assert!(cmd(&["cache", "status"]).status.success());
    assert!(!root.exists());
    assert!(cmd(&["cache", "prune", "--all-evictable", "--dry-run"])
        .status
        .success());
    assert!(!root.exists());
    let mut c = Cache::open(&root, DEFAULT_LIMIT, true).unwrap();
    let mut m = manifest();
    put(&mut c, &mut m, b"keep", 0, Retention::Keep);
    put(&mut c, &mut m, b"discard", 0, Retention::Evict);
    drop(c);
    let r = cmd(&["cache", "prune", "--all-evictable", "--dry-run"]);
    assert!(r.status.success());
    let value: serde_json::Value = serde_json::from_slice(&r.stdout).unwrap();
    assert_eq!(value["freed_bytes"], 7);
    assert!(cmd(&["cache", "prune", "--all-evictable"]).status.success());
    let r = cmd(&["cache", "status"]);
    let value: serde_json::Value = serde_json::from_slice(&r.stdout).unwrap();
    assert_eq!(value["large_data_bytes"], 4);
    assert!(cmd(&["cache", "purge-data"]).status.success());
    let r = cmd(&["cache", "status"]);
    let value: serde_json::Value = serde_json::from_slice(&r.stdout).unwrap();
    assert_eq!(value["large_data_bytes"], 4);
    assert!(cmd(&["cache", "purge-data", "--yes"]).status.success());
    let r = cmd(&["manifest", "PXD000001"]);
    assert!(r.status.success());
    let r = cmd(&["cache", "status"]);
    let value: serde_json::Value = serde_json::from_slice(&r.stdout).unwrap();
    assert_eq!(value["large_data_bytes"], 0);
}

#[test]
fn remote_identity_ignores_display_metadata_and_reference_order() {
    let mut f = file("https://example.org/a", b"abc");
    f.checksums.clear();
    f.references.push("https://example.org/b".into());
    let key = f.object_key();
    f.references.reverse();
    f.category = Some("OTHER".into());
    f.filename = "new-display-name.pin".into();
    assert_eq!(key, f.object_key());
    f.size_bytes = Some(4);
    assert_ne!(key, f.object_key());
}
#[test]
fn cached_inputs_still_count_toward_working_space_and_staging_is_reported() {
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 100);
    let mut m = manifest();
    put(&mut c, &mut m, b"12345678", 0, Retention::Evict);
    let f = file("https://example.org/a", b"12345678");
    let mut b = budgets();
    b.max_download = 0;
    b.max_working_space = Some(9);
    assert!(download::plan(&c, &m, std::slice::from_ref(&f), &b, false, 1, 4).is_err());
    b.max_working_space = Some(12);
    let p = download::plan(&c, &m, &[f], &b, true, 1, 4).unwrap();
    assert_eq!(p.download_bytes, 0);
    assert_eq!(p.temporary_workspace_bytes, 4);
    assert_eq!(p.peak_working_bytes, 12);
    assert_eq!(p.expected_final_large_data_bytes, 0);
}

#[test]
fn failed_preparation_preserves_recipe_without_claiming_a_prepared_pin() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("invalid.pin");
    fs::write(&input, b"invalid PIN\n").unwrap();
    let mut c = cache(&tmp, 100);
    let mut m = manifest();
    m.inventory
        .push(file("https://example.org/source", b"source"));
    assert!(pride::prepare::import(
        &mut c,
        &mut m,
        &input,
        recipe("remote-id", "bad-export"),
        &budgets(),
        Retention::KeepIfPinned
    )
    .is_err());
    let m = c.load_manifest(&pxd()).unwrap();
    assert!(m.prepared_pins.is_empty());
    assert_eq!(m.preparation_attempts[0].state, "failed");
    assert_eq!(m.preparation_attempts[0].steps[0].inputs, vec!["remote-id"]);
    c.prune(true, false, false).unwrap();
    assert_eq!(
        c.load_manifest(&pxd()).unwrap().preparation_attempts.len(),
        1
    );
}

#[test]
fn resumed_operation_reserves_only_remaining_bytes_under_a_full_cache_ceiling() {
    let healthy = Arc::new(AtomicBool::new(false));
    let flag = healthy.clone();
    let server = Server::new(move |r| {
        if flag.load(Ordering::Relaxed) {
            assert!(r.to_ascii_lowercase().contains("range: bytes=3-"));
            Reply {
                status: 206,
                headers: vec![("Content-Range".into(), "bytes 3-5/6".into())],
                ..Reply::ok("def")
            }
        } else if r.to_ascii_lowercase().contains("range:") {
            Reply {
                status: 503,
                ..Reply::ok("busy")
            }
        } else {
            Reply {
                length: Some(6),
                ..Reply::ok("abc")
            }
        }
    });
    let tmp = TempDir::new().unwrap();
    let mut c = cache(&tmp, 6);
    let mut m = manifest();
    let f = file(&server.url, b"abcdef");
    assert!(downloader()
        .fetch(&mut c, &mut m, &f, &BTreeSet::new())
        .is_err());
    assert_eq!(c.status().unwrap().large_data_bytes, 3);
    let root = c.root.clone();
    drop(c);
    let mut c = Cache::open(&root, 6, true).unwrap();
    let mut b = budgets();
    b.max_download = 3;
    let p = download::plan(&c, &m, std::slice::from_ref(&f), &b, false, 1, 0).unwrap();
    assert_eq!(p.download_bytes, 3);
    assert_eq!(p.expected_final_large_data_bytes, 6);
    healthy.store(true, Ordering::Relaxed);
    let mut dl = Downloader::new(&b, Arc::new(AtomicBool::new(false))).unwrap();
    dl.fetch(&mut c, &mut m, &f, &BTreeSet::new()).unwrap();
    assert_eq!(dl.remaining_download, 0);
    assert_eq!(c.status().unwrap().large_data_bytes, 6);
    assert_eq!(c.index.objects[&f.object_key()].state, State::Verified);
}
