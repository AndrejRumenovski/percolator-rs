//! The only module aware of PRIDE v3 wire formats.
use super::*;
use reqwest::{
    blocking::{Client, Response},
    Url,
};
use serde::Deserialize;
use serde_json::Value;
use std::{collections::HashSet, io::Read, time::Duration};

const METADATA_LIMIT: u64 = 32 * 1024 * 1024;
#[derive(Clone)]
pub struct PrideClient {
    http: Client,
    base: String,
}
impl PrideClient {
    pub fn new() -> Result<Self> {
        Self::with_base(API_BASE)
    }
    /// Alternate origins are for controlled mirrors/tests; the CLI always uses official PRIDE.
    pub fn with_base(base: &str) -> Result<Self> {
        let url = Url::parse(base)?;
        if !matches!(url.scheme(), "http" | "https") {
            return Err("PRIDE API requires HTTP(S)".into());
        }
        Ok(Self {
            http: Client::builder()
                .connect_timeout(Duration::from_secs(15))
                .timeout(Duration::from_secs(90))
                .user_agent(concat!(
                    "percolator-rs/",
                    env!("CARGO_PKG_VERSION"),
                    " PRIDE-client"
                ))
                .redirect(reqwest::redirect::Policy::limited(5))
                .build()?,
            base: base.trim_end_matches('/').into(),
        })
    }
    fn get(&self, path: &str, query: &[(&str, String)]) -> Result<Response> {
        let mut last = String::new();
        for attempt in 0..3 {
            match self
                .http
                .get(format!("{}{path}", self.base))
                .query(query)
                .send()
            {
                Ok(r) if r.status().is_success() => return Ok(r),
                Ok(r) => {
                    let status = r.status();
                    last = format!(
                        "PRIDE {path}: HTTP {status}{}",
                        if status.as_u16() == 404 {
                            " (project/file unavailable; check accession and public status)"
                        } else {
                            ""
                        }
                    );
                    if !status.is_server_error() && status.as_u16() != 429 {
                        return Err(last.into());
                    }
                }
                Err(e) => last = format!("PRIDE {path}: {e}"),
            }
            if attempt < 2 {
                std::thread::sleep(Duration::from_millis(250 << attempt));
            }
        }
        Err(last.into())
    }
    fn text(&self, path: &str, query: &[(&str, String)]) -> Result<String> {
        let mut body = String::new();
        self.get(path, query)?
            .take(METADATA_LIMIT + 1)
            .read_to_string(&mut body)?;
        if body.len() as u64 > METADATA_LIMIT {
            return Err("PRIDE metadata response exceeds 32 MiB safety limit".into());
        }
        Ok(body)
    }
    fn json<T: serde::de::DeserializeOwned>(
        &self,
        path: &str,
        query: &[(&str, String)],
    ) -> Result<T> {
        serde_json::from_str(&self.text(path, query)?)
            .map_err(|e| format!("PRIDE {path}: invalid v3 response: {e}").into())
    }
    pub fn project(&self, accession: &Pxd) -> Result<Project> {
        let status = self.text(&format!("/status/{accession}"), &[])?;
        let status = status.trim().trim_matches('"');
        if status != "PUBLIC" {
            return Err(format!(
                "{accession} is not accessible as a public PRIDE project (status {status:?})"
            )
            .into());
        }
        let mut p = parse_project(self.json(&format!("/projects/{accession}"), &[])?)?;
        if &p.accession != accession {
            return Err("PRIDE returned the wrong project accession".into());
        }
        p.status = Some(status.into());
        Ok(p)
    }
    pub fn files(&self, accession: &Pxd) -> Result<Vec<RemoteFile>> {
        let count_path = format!("/projects/{accession}/files/count");
        let expected: usize = self.json(&count_path, &[])?;
        if expected > 1_000_000 {
            return Err(
                "PRIDE inventory exceeds one million records; refusing unbounded metadata".into(),
            );
        }
        let mut files = Vec::new();
        let mut seen = HashSet::new();
        let mut page = 0;
        let mut metadata_bytes = 0u64;
        while files.len() < expected {
            let raw: Vec<WireFile> = self.json(
                &format!("/projects/{accession}/files"),
                &[("pageSize", "100".into()), ("page", page.to_string())],
            )?;
            if raw.is_empty() {
                return Err(format!(
                    "incomplete PRIDE inventory: {} of {expected} records",
                    files.len()
                )
                .into());
            }
            for f in raw {
                let f = f.into_domain(&self.base)?;
                if !seen.insert(f.id.clone()) {
                    return Err(
                        "PRIDE repeated a file/page; refusing an incomplete or unstable inventory"
                            .into(),
                    );
                }
                metadata_bytes = total([metadata_bytes, serde_json::to_vec(&f)?.len() as u64])?;
                if metadata_bytes > METADATA_LIMIT {
                    return Err(
                        "complete PRIDE inventory exceeds 32 MiB metadata safety limit".into(),
                    );
                }
                files.push(f);
            }
            page += 1;
        }
        let after: usize = self.json(&count_path, &[])?;
        if files.len() != expected || after != expected {
            return Err(
                "PRIDE file inventory changed while paging; retry metadata retrieval".into(),
            );
        }
        Ok(files)
    }
    pub fn manifest(&self, accession: &Pxd) -> Result<Manifest> {
        let project = self.project(accession)?;
        let mut m = Manifest::new(project, self.files(accession)?, self.base.clone());
        let endpoint = format!("/files/checksum/{accession}");
        // A missing checksum table is a legitimate absence; other HTTP/parsing errors
        // are surfaced, never silently converted into 'no checksum'.
        let table = match self.text(&endpoint, &[]) {
            Ok(t) => Some(t),
            Err(e) if e.to_string().contains("HTTP 404") => None,
            Err(e) => return Err(e),
        };
        if let Some(table) = table {
            let paths: Value = self.json(&format!("/projects/files-path/{accession}"), &[])?;
            let root = paths.get("ftp").and_then(Value::as_str);
            merge_checksums(&mut m, &table, &format!("{}{endpoint}", self.base), root)?;
        } else {
            m.inventory_notes
                .push("Repository checksum table is unavailable (404).".into());
        }
        m.inventory_notes.push("Inventory covers every indexed API record plus resolvable checksum-table entries; PRIDE may hold additional unindexed archive files.".into());
        Ok(m)
    }
    pub fn search(
        &self,
        keyword: &str,
        filter: Option<&str>,
        page: u32,
        page_size: u32,
    ) -> Result<Vec<Project>> {
        if !(1..=100).contains(&page_size) {
            return Err("page size must be 1..100".into());
        }
        let mut q = vec![
            ("keyword", keyword.into()),
            ("page", page.to_string()),
            ("pageSize", page_size.to_string()),
        ];
        if let Some(f) = filter {
            q.push(("filter", f.into()));
        }
        let rows: Vec<Value> = self.json("/search/projects", &q)?;
        rows.into_iter().map(parse_project).collect()
    }
}
impl Default for PrideClient {
    fn default() -> Self {
        Self::new().expect("default HTTP client")
    }
}

fn opt_string(v: &Value, key: &str) -> Result<Option<String>> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.clone())),
        _ => Err(format!("PRIDE {key} must be a string or null").into()),
    }
}
fn terms(v: &Value, key: &str) -> Result<Option<Vec<Term>>> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(a)) => Ok(Some(
            a.iter()
                .map(|v| match v {
                    Value::String(s) => Ok(Term {
                        name: Some(s.clone()),
                        ..Term::default()
                    }),
                    Value::Object(_) => serde_json::from_value(v.clone()).map_err(Into::into),
                    _ => Err(format!("invalid PRIDE {key} term").into()),
                })
                .collect::<Result<_>>()?,
        )),
        _ => Err(format!("PRIDE {key} must be an array or null").into()),
    }
}
pub fn parse_project(v: Value) -> Result<Project> {
    let accession = opt_string(&v, "accession")?
        .ok_or("project accession missing")?
        .parse()?;
    let publications = match v.get("references") {
        None | Some(Value::Null) => None,
        Some(Value::Array(a)) => Some(
            a.iter()
                .map(|r| -> Result<Publication> {
                    if let Some(s) = r.as_str() {
                        return Ok(Publication {
                            citation: Some(s.into()),
                            ..Publication::default()
                        });
                    }
                    Ok(Publication {
                        citation: opt_string(r, "referenceLine")?,
                        doi: opt_string(r, "doi")?,
                        pubmed_id: r.get("pubmedID").filter(|v| !v.is_null()).map(|v| {
                            v.as_str()
                                .map(str::to_owned)
                                .unwrap_or_else(|| v.to_string())
                        }),
                    })
                })
                .collect::<Result<_>>()?,
        ),
        _ => return Err("invalid PRIDE references".into()),
    };
    Ok(Project {
        accession,
        title: opt_string(&v, "title")?,
        description: opt_string(&v, "projectDescription")?,
        organisms: terms(&v, "organisms")?,
        tissues: terms(&v, "organismParts")?.or(terms(&v, "organismsPart")?),
        instruments: terms(&v, "instruments")?,
        experiment_types: terms(&v, "experimentTypes")?,
        modifications: terms(&v, "identifiedPTMStrings")?.or(terms(&v, "identifiedPtmStrings")?),
        publications,
        submission_date: opt_string(&v, "submissionDate")?,
        publication_date: opt_string(&v, "publicationDate")?,
        submission_type: opt_string(&v, "submissionType")?,
        doi: opt_string(&v, "doi")?,
        sample_processing_protocol: opt_string(&v, "sampleProcessingProtocol")?,
        data_processing_protocol: opt_string(&v, "dataProcessingProtocol")?,
        status: None,
    })
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct WireFile {
    accession: String,
    file_name: String,
    file_category: Option<Term>,
    file_extension: Option<String>,
    file_size_bytes: Option<u64>,
    public_file_locations: Option<Vec<Term>>,
    checksum: Option<String>,
    analysis_accessions: Option<Vec<String>>,
    additional_attributes: Option<Vec<Term>>,
}
impl WireFile {
    fn into_domain(self, base: &str) -> Result<RemoteFile> {
        if self.accession.is_empty() || self.file_name.is_empty() {
            return Err("PRIDE file identity/name missing".into());
        }
        let raw = self.checksum.filter(|s| !s.is_empty());
        let mut checksums = Vec::new();
        let mut untyped_checksum = raw.clone();
        // PRIDE's current checksum generator is SHA-1. Do not interpret ambiguous
        // shorter values as MD5: only the separately labelled legacy table is MD5.
        if let Some(s) = &raw {
            if s.len() == 40 {
                checksums.push(Checksum::new(
                    Algorithm::Sha1,
                    s,
                    &format!("{base}/files/{}#checksum", self.accession),
                )?);
                untyped_checksum = None;
            }
        }
        Ok(RemoteFile {
            id: self.accession,
            filename: self.file_name,
            category: self.file_category.and_then(|t| t.value),
            format: self.file_extension,
            size_bytes: self.file_size_bytes,
            checksum_table_size: None,
            references: self
                .public_file_locations
                .unwrap_or_default()
                .into_iter()
                .filter_map(|t| t.value)
                .collect(),
            checksums,
            untyped_checksum,
            analysis_accessions: self.analysis_accessions,
            run_metadata: self.additional_attributes,
            inventory_source: "archive_v3_files".into(),
        })
    }
}
pub fn parse_files(v: Value) -> Result<Vec<RemoteFile>> {
    serde_json::from_value::<Vec<WireFile>>(v)?
        .into_iter()
        .map(|f| f.into_domain(API_BASE))
        .collect()
}

/// Preserve both size claims and refuse ambiguous name matches. Supplemental entries
/// use the official project FTP root, never a date inferred from project metadata.
pub fn merge_checksums(
    m: &mut Manifest,
    text: &str,
    authority: &str,
    root: Option<&str>,
) -> Result<()> {
    if text.trim().is_empty() {
        m.inventory_notes
            .push("Empty repository checksum table.".into());
        return Ok(());
    }
    let mut lines = text.lines();
    let header = lines.next().unwrap_or("").trim_end_matches('\r');
    if header != "File-Name\tFile-MD5Checksum\tFile-Size" {
        return Err(format!("unsupported PRIDE checksum-table header {header:?}").into());
    }
    let mut seen = HashSet::new();
    let mut metadata_bytes = serde_json::to_vec(&m.inventory)?.len() as u64;
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let cols: Vec<_> = line.trim_end_matches('\r').split('\t').collect();
        if cols.len() != 3 || !seen.insert(cols[0]) {
            return Err("malformed or duplicate PRIDE checksum-table row".into());
        }
        let size: u64 = cols[2].parse()?;
        let checksum = if cols[1].is_empty() || cols[1] == "null" {
            None
        } else {
            Some(Checksum::new(Algorithm::Md5, cols[1], authority)?)
        };
        let matches: Vec<_> = m
            .inventory
            .iter()
            .enumerate()
            .filter(|(_, f)| f.filename == cols[0])
            .map(|(i, _)| i)
            .collect();
        if matches.len() > 1 {
            return Err(format!("ambiguous checksum filename {:?}", cols[0]).into());
        }
        if let Some(&i) = matches.first() {
            let f = &mut m.inventory[i];
            f.checksum_table_size = Some(size);
            if let Some(c) = checksum {
                f.checksums.push(c);
            }
        } else if let Some(root) = root {
            // This is a remote relative name only; never used as a local path.
            if !safe_remote_relative(cols[0]) {
                return Err("unsafe checksum-table filename".into());
            }
            let mut u = Url::parse(root)?;
            {
                let mut p = u
                    .path_segments_mut()
                    .map_err(|_| "invalid PRIDE FTP root")?;
                p.pop_if_empty();
                for part in cols[0].split('/') {
                    p.push(part);
                }
            }
            metadata_bytes = total([
                metadata_bytes,
                u.as_str().len() as u64,
                cols[0].len() as u64 * 2,
                1024,
            ])?;
            if metadata_bytes > METADATA_LIMIT {
                return Err("supplemental PRIDE inventory exceeds metadata safety limit".into());
            }
            m.inventory.push(RemoteFile {
                id: format!("checksum-table:{}", cols[0]),
                filename: cols[0].into(),
                category: None,
                format: None,
                size_bytes: None,
                checksum_table_size: Some(size),
                references: vec![u.to_string()],
                checksums: checksum.into_iter().collect(),
                untyped_checksum: None,
                analysis_accessions: None,
                run_metadata: None,
                inventory_source: "archive_v3_checksum_table".into(),
            });
        } else {
            m.inventory_notes.push(format!(
                "Checksum entry {:?} has no resolvable project root",
                cols[0]
            ));
        }
    }
    Ok(())
}
pub fn safe_remote_relative(name: &str) -> bool {
    !name.is_empty()
        && !name.contains(['\\', ':'])
        && !name.chars().any(char::is_control)
        && name
            .split('/')
            .all(|s| !s.is_empty() && s != "." && s != "..")
}
/// Only the documented PRIDE FTP -> HTTPS mapping is implicit. Other HTTPS
/// references remain usable; Aspera/FTP on other hosts requires external tooling.
pub fn download_url(file: &RemoteFile) -> Result<Url> {
    for reference in &file.references {
        if let Ok(mut u) = Url::parse(reference) {
            if u.scheme() == "ftp"
                && u.host_str() == Some("ftp.pride.ebi.ac.uk")
                && u.port().is_none()
            {
                u.set_scheme("https")
                    .map_err(|_| "cannot resolve PRIDE HTTPS URL")?;
            }
            if u.scheme() == "https" && u.username().is_empty() && u.password().is_none() {
                return Ok(u);
            }
            // Deliberately narrow support for local offline HTTP fixtures.
            if u.scheme() == "http" && matches!(u.host_str(), Some("127.0.0.1" | "[::1]")) {
                return Ok(u);
            }
        }
    }
    Err(format!("{}: no supported HTTPS download reference", file.filename).into())
}
