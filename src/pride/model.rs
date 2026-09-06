use super::{now, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, fmt, str::FromStr};

pub const MANIFEST_SCHEMA_VERSION: u32 = 1;
pub const API_BASE: &str = "https://www.ebi.ac.uk/pride/ws/archive/v3";

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct Pxd(String);
impl FromStr for Pxd {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        if s.len() == 9 && s.starts_with("PXD") && s[3..].bytes().all(|b| b.is_ascii_digit()) {
            Ok(Self(s.to_owned()))
        } else {
            Err(format!("invalid PRIDE accession {s:?}: expected PXD followed by exactly six digits (e.g. PXD000001)"))
        }
    }
}
impl TryFrom<String> for Pxd {
    type Error = String;
    fn try_from(s: String) -> std::result::Result<Self, String> {
        s.parse()
    }
}
impl From<Pxd> for String {
    fn from(p: Pxd) -> String {
        p.0
    }
}
impl fmt::Display for Pxd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Term {
    pub accession: Option<String>,
    pub name: Option<String>,
    pub value: Option<String>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Publication {
    pub citation: Option<String>,
    pub doi: Option<String>,
    pub pubmed_id: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Project {
    pub accession: Pxd,
    pub title: Option<String>,
    pub description: Option<String>,
    pub organisms: Option<Vec<Term>>,
    pub tissues: Option<Vec<Term>>,
    pub instruments: Option<Vec<Term>>,
    pub experiment_types: Option<Vec<Term>>,
    pub modifications: Option<Vec<Term>>,
    pub publications: Option<Vec<Publication>>,
    pub submission_date: Option<String>,
    pub publication_date: Option<String>,
    pub submission_type: Option<String>,
    pub doi: Option<String>,
    pub sample_processing_protocol: Option<String>,
    pub data_processing_protocol: Option<String>,
    pub status: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    Md5,
    Sha1,
    Sha256,
}
impl Algorithm {
    pub fn label(self) -> &'static str {
        match self {
            Self::Md5 => "md5",
            Self::Sha1 => "sha1",
            Self::Sha256 => "sha256",
        }
    }
    pub fn width(self) -> usize {
        match self {
            Self::Md5 => 32,
            Self::Sha1 => 40,
            Self::Sha256 => 64,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Checksum {
    pub algorithm: Algorithm,
    pub value: String,
    pub reported_value: String,
    pub authority: String,
}
impl Checksum {
    pub fn new(algorithm: Algorithm, value: &str, authority: &str) -> Result<Self> {
        if value.is_empty()
            || value.len() != algorithm.width()
            || !value.bytes().all(|c| c.is_ascii_hexdigit())
        {
            return Err(format!(
                "invalid authoritative {} checksum {value:?}",
                algorithm.label()
            )
            .into());
        }
        Ok(Self {
            algorithm,
            value: format!(
                "{:0>width$}",
                value.to_ascii_lowercase(),
                width = algorithm.width()
            ),
            reported_value: value.into(),
            authority: authority.into(),
        })
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Compatibility {
    DirectlyCompatible,
    PotentiallyConvertible,
    RawRequiresSearch,
    UnrelatedUnknown,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RemoteFile {
    pub id: String,
    pub filename: String,
    pub category: Option<String>,
    pub format: Option<String>,
    pub size_bytes: Option<u64>,
    pub checksum_table_size: Option<u64>,
    pub references: Vec<String>,
    pub checksums: Vec<Checksum>,
    pub untyped_checksum: Option<String>,
    pub analysis_accessions: Option<Vec<String>>,
    pub run_metadata: Option<Vec<Term>>,
    pub inventory_source: String,
}
impl RemoteFile {
    pub fn format_name(&self) -> String {
        let name = self.filename.to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        if name.ends_with(".pep.xml") {
            "pepxml".into()
        } else if name.ends_with(".mzid") || name.ends_with(".mzidentml") {
            "mzidentml".into()
        } else if name.ends_with("mztab.txt") {
            "mztab".into()
        } else {
            name.rsplit('.').next().unwrap_or("").into()
        }
    }
    pub fn native_pin(&self) -> bool {
        self.filename.to_ascii_lowercase().ends_with(".pin")
    }
    pub fn compatibility(&self) -> Compatibility {
        match self.format_name().as_str() {
            "pin" | "mzidentml" | "mztab" | "pepxml" | "sqt" | "dat" => {
                Compatibility::PotentiallyConvertible
            }
            "raw" | "mzml" | "mzxml" | "mgf" | "wiff" | "d" | "tdf" | "baf" => {
                Compatibility::RawRequiresSearch
            }
            _ if self.category.as_deref() == Some("RAW") => Compatibility::RawRequiresSearch,
            _ if self.category.as_deref() == Some("SEARCH") => {
                Compatibility::PotentiallyConvertible
            }
            _ => Compatibility::UnrelatedUnknown,
        }
    }
    pub fn preparation(&self) -> &'static str {
        if self.native_pin() {
            return "Native PIN candidate: existing parser must validate content before direct compatibility is confirmed; the concatenated target/decoy search contract still applies";
        }
        match self.compatibility() {
            Compatibility::DirectlyCompatible => "PIN candidate: validate with the existing parser; requires the documented concatenated target/decoy search contract",
            Compatibility::PotentiallyConvertible => "External preparation required: export full target/decoy candidates and numeric features to PIN (or decompress PIN); no automatic converter installed",
            Compatibility::RawRequiresSearch => "External conversion if needed, database search with an explicit target/decoy database and parameters, then PIN export",
            Compatibility::UnrelatedUnknown => "Auxiliary/unknown file; inspect metadata before selecting a workflow",
        }
    }
    pub fn matches(&self, filter: &str) -> bool {
        let f = filter.to_ascii_lowercase();
        match f.as_str() {
            "processed" => {
                matches!(self.category.as_deref(), Some("RESULT" | "SEARCH"))
                    || matches!(
                        self.compatibility(),
                        Compatibility::DirectlyCompatible | Compatibility::PotentiallyConvertible
                    )
            }
            "search-engine-output" => {
                self.category.as_deref() == Some("SEARCH")
                    || matches!(self.format_name().as_str(), "pepxml" | "sqt" | "dat")
            }
            "mzid" => self.format_name() == "mzidentml",
            _ => {
                self.category
                    .as_ref()
                    .is_some_and(|x| x.eq_ignore_ascii_case(&f))
                    || self.format_name() == f
                    || self
                        .format
                        .as_ref()
                        .is_some_and(|x| x.eq_ignore_ascii_case(&f))
            }
        }
    }
    pub fn size(&self) -> Result<u64> {
        if let (Some(a), Some(b)) = (self.size_bytes, self.checksum_table_size) {
            if a != b {
                return Err(format!("{} has contradictory API/table sizes ({a}/{b}); refresh or select another file", self.filename).into());
            }
        }
        self.size_bytes.or(self.checksum_table_size).ok_or_else(|| {
            format!(
                "{}: remote size missing; cannot safely budget this download",
                self.filename
            )
            .into()
        })
    }
    pub fn object_key(&self) -> String {
        // An authoritative identity enables sharing before downloading. Otherwise bind the
        // remote record, all references, size, and checksum evidence to a local namespace.
        if let Some(c) = self.checksums.iter().max_by_key(|c| c.algorithm.width()) {
            format!("{}-{}", c.algorithm.label(), c.value)
        } else {
            let mut references = self.references.clone();
            references.sort();
            references.dedup();
            let identity = (
                &self.id,
                references,
                self.size_bytes,
                self.checksum_table_size,
                &self.untyped_checksum,
            );
            format!(
                "remote-{:x}",
                Sha256::digest(
                    serde_json::to_vec(&identity).expect("serializable remote identity")
                )
            )
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum State {
    Remote,
    Downloading,
    Partial,
    Verified,
    DownloadedUnverified,
    Prepared,
    Corrupt,
    Evicted,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Retention {
    Keep,
    Evict,
    KeepIfPinned,
    UntilResultVerified,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Verification {
    pub algorithm: Algorithm,
    pub expected: String,
    pub actual: String,
    pub verified: bool,
    pub authority: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalFile {
    pub object_key: String,
    pub local_relative_path: String,
    pub state: State,
    pub local_sha256: Option<String>,
    pub verification: Vec<Verification>,
    pub verification_unavailable: bool,
    pub pin_validated: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lineage {
    pub id: String,
    pub inputs: Vec<String>,
    pub output_sha256: Option<String>,
    pub kind: String,
    pub tool: String,
    pub tool_version: Option<String>,
    pub parameters: Vec<String>,
    pub protein_database: Option<String>,
    pub database_sha256: Option<String>,
    pub decoy_generation: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Experiment {
    pub id: String,
    pub input_ids: Vec<String>,
    pub state: String,
    pub error: Option<String>,
    pub started_unix_seconds: u64,
    pub completed_unix_seconds: Option<u64>,
    pub percolator_rs_version: String,
    pub percolator_rs_commit: Option<String>,
    pub executable_sha256: String,
    pub parameters: Vec<String>,
    pub ephemeral: bool,
    pub pin_retention: Retention,
    pub result_hashes: BTreeMap<String, String>,
    pub lineage: Vec<Lineage>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreparationAttempt {
    pub steps: Vec<Lineage>,
    pub state: String,
    pub error: Option<String>,
    pub started_unix_seconds: u64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreparedPin {
    pub id: String,
    pub object_key: String,
    pub sha256: String,
    pub bytes: u64,
    pub lineage_id: String,
    pub retention: Retention,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub schema_version: u32,
    pub accession: Pxd,
    pub project: Project,
    pub api_base: String,
    pub retrieved_unix_seconds: u64,
    pub indexed_file_count: usize,
    pub inventory: Vec<RemoteFile>,
    pub remote_history: Vec<RemoteFile>,
    pub prepared_pins: BTreeMap<String, PreparedPin>,
    #[serde(default)]
    pub preparation_attempts: Vec<PreparationAttempt>,
    pub inventory_notes: Vec<String>,
    pub selected_files: Vec<String>,
    pub local_files: BTreeMap<String, LocalFile>,
    pub experiments: Vec<Experiment>,
    pub lineage: Vec<Lineage>,
}
impl Manifest {
    pub fn new(project: Project, inventory: Vec<RemoteFile>, api_base: String) -> Self {
        Self {
            schema_version: MANIFEST_SCHEMA_VERSION,
            accession: project.accession.clone(),
            project,
            api_base,
            retrieved_unix_seconds: now(),
            indexed_file_count: inventory.len(),
            inventory,
            remote_history: vec![],
            prepared_pins: BTreeMap::new(),
            preparation_attempts: vec![],
            inventory_notes: vec![],
            selected_files: vec![],
            local_files: BTreeMap::new(),
            experiments: vec![],
            lineage: vec![],
        }
    }
    pub fn compatibility(&self, f: &RemoteFile) -> Compatibility {
        if f.native_pin()
            && self.local_files.get(&f.id).is_some_and(|l| {
                l.pin_validated && l.object_key == f.object_key() && l.state != State::Corrupt
            })
        {
            Compatibility::DirectlyCompatible
        } else {
            f.compatibility()
        }
    }
    pub fn validate(&self) -> Result<()> {
        if self.schema_version != MANIFEST_SCHEMA_VERSION
            || self.accession != self.project.accession
            || self.indexed_file_count > self.inventory.len()
        {
            return Err(
                "unsupported or inconsistent PRIDE manifest; preserved without modification".into(),
            );
        }
        let mut ids = std::collections::HashSet::new();
        for f in &self.inventory {
            if !ids.insert(&f.id) {
                return Err("duplicate file identities in manifest".into());
            }
        }
        Ok(())
    }
}
