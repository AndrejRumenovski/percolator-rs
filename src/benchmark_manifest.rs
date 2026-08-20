//! Dataset metadata used by benchmark tooling.
//!
//! This module deliberately only describes and validates benchmark inputs. It
//! does not select datasets, discover files, download data, or run Percolator.

use serde::Deserialize;
use std::fmt;
use std::path::Path;

/// Top-level dataset registry stored in `bench/datasets.toml`.
#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DatasetRegistry {
    pub version: u32,
    #[serde(rename = "datasets")]
    pub datasets: Vec<Dataset>,
}

/// Metadata and the local PIN location for one benchmark dataset.
#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Dataset {
    pub id: String,
    #[serde(default)]
    pub pride_accession: Option<String>,
    pub source: String,
    pub organism: String,
    pub experiment_type: String,
    #[serde(default)]
    pub instrument: Option<String>,
    pub search_engine: String,
    pub pin_path: String,
    #[serde(default)]
    pub file_count: Option<usize>,
    #[serde(default)]
    pub approximate_input_size: Option<String>,
    pub protein_level_evaluation: bool,
    pub notes: String,
    #[serde(default)]
    pub preparation: Option<String>,
    /// C++ Percolator's explicit target/decoy input interpretation, when the
    /// dataset requires one. This is deliberately dataset metadata rather
    /// than a runner default because it changes the statistical methodology.
    #[serde(default)]
    pub reference_search_input: Option<SearchInput>,
}

/// Supported values for C++ Percolator's `--search-input` option.
#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum SearchInput {
    Concatenated,
    Separate,
}

impl SearchInput {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Concatenated => "concatenated",
            Self::Separate => "separate",
        }
    }
}

/// An actionable problem found while loading a benchmark manifest.
#[derive(Debug, PartialEq, Eq)]
pub enum ManifestError {
    Parse(String),
    Validation(Vec<String>),
}

impl fmt::Display for ManifestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(message) => write!(f, "invalid dataset manifest: {message}"),
            Self::Validation(errors) => {
                writeln!(f, "invalid dataset manifest:")?;
                for error in errors {
                    writeln!(f, "  - {error}")?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for ManifestError {}

impl DatasetRegistry {
    /// Parse and validate TOML manifest contents.
    pub fn from_toml(input: &str) -> Result<Self, ManifestError> {
        let registry: Self =
            toml::from_str(input).map_err(|error| ManifestError::Parse(error.to_string()))?;
        registry.validate()?;
        Ok(registry)
    }

    /// Read and validate a registry from disk.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, ManifestError> {
        let path = path.as_ref();
        let input = std::fs::read_to_string(path).map_err(|error| {
            ManifestError::Parse(format!("could not read {}: {error}", path.display()))
        })?;
        Self::from_toml(&input)
    }

    /// Validate cross-dataset and individual metadata constraints.
    pub fn validate(&self) -> Result<(), ManifestError> {
        let mut errors = Vec::new();
        if self.version != 1 {
            errors.push(format!("unsupported version {}; expected 1", self.version));
        }
        if self.datasets.is_empty() {
            errors.push("at least one [[datasets]] entry is required".to_owned());
        }

        let mut ids = std::collections::BTreeSet::new();
        for dataset in &self.datasets {
            let label = format!("dataset {:?}", dataset.id);
            if !is_dataset_id(&dataset.id) {
                errors.push(format!("{label}: id must use ASCII letters, digits, '-' or '_' and start with a letter or digit"));
            }
            if !ids.insert(&dataset.id) {
                errors.push(format!("{label}: duplicate id"));
            }
            require_text(&mut errors, &label, "source", &dataset.source);
            require_text(&mut errors, &label, "organism", &dataset.organism);
            require_text(
                &mut errors,
                &label,
                "experiment_type",
                &dataset.experiment_type,
            );
            require_text(&mut errors, &label, "search_engine", &dataset.search_engine);
            require_text(&mut errors, &label, "pin_path", &dataset.pin_path);
            require_text(&mut errors, &label, "notes", &dataset.notes);

            if !dataset.pin_path.ends_with(".pin") {
                errors.push(format!(
                    "{label}: pin_path must name a .pin file or glob ending in .pin"
                ));
            }
            if let Err(error) = validate_environment_templates(&dataset.pin_path) {
                errors.push(format!("{label}: pin_path {error}"));
            }
            if dataset.file_count == Some(0) {
                errors.push(format!(
                    "{label}: file_count must be greater than zero when provided"
                ));
            }
            optional_text(
                &mut errors,
                &label,
                "pride_accession",
                &dataset.pride_accession,
            );
            optional_text(&mut errors, &label, "instrument", &dataset.instrument);
            optional_text(
                &mut errors,
                &label,
                "approximate_input_size",
                &dataset.approximate_input_size,
            );
            optional_text(&mut errors, &label, "preparation", &dataset.preparation);
        }

        if errors.is_empty() {
            Ok(())
        } else {
            Err(ManifestError::Validation(errors))
        }
    }
}

fn is_dataset_id(id: &str) -> bool {
    let mut chars = id.bytes();
    matches!(chars.next(), Some(byte) if byte.is_ascii_alphanumeric())
        && chars.all(|byte| byte.is_ascii_alphanumeric() || byte == b'-' || byte == b'_')
}

fn require_text(errors: &mut Vec<String>, label: &str, field: &str, value: &str) {
    if value.trim().is_empty() {
        errors.push(format!("{label}: {field} must not be empty"));
    }
}

fn optional_text(errors: &mut Vec<String>, label: &str, field: &str, value: &Option<String>) {
    if value.as_deref().is_some_and(|text| text.trim().is_empty()) {
        errors.push(format!("{label}: {field} must not be empty when provided"));
    }
}

/// Ensure `${NAME}` templates are complete and use portable environment names.
/// Expansion is deliberately deferred to a future benchmark runner.
fn validate_environment_templates(value: &str) -> Result<(), &'static str> {
    let mut remainder = value;
    while let Some(start) = remainder.find("${") {
        let after_start = &remainder[start + 2..];
        let Some(end) = after_start.find('}') else {
            return Err("contains an unterminated ${...} environment template");
        };
        let name = &after_start[..end];
        if name.is_empty()
            || !name
                .bytes()
                .all(|byte| byte.is_ascii_uppercase() || byte.is_ascii_digit() || byte == b'_')
        {
            return Err("contains an invalid ${...} environment template; use an uppercase name such as ${PERCOLATOR_BENCH_DATA}");
        }
        remainder = &after_start[end + 1..];
    }
    Ok(())
}

/// Expand `${UPPERCASE_ENV}` path templates using the process environment.
///
/// The manifest validates template syntax on load; this function reports a
/// useful runtime error if a required machine-specific location is unset.
pub fn expand_environment_templates(value: &str) -> Result<String, String> {
    validate_environment_templates(value).map_err(str::to_owned)?;
    let mut expanded = String::with_capacity(value.len());
    let mut remainder = value;
    while let Some(start) = remainder.find("${") {
        expanded.push_str(&remainder[..start]);
        let after_start = &remainder[start + 2..];
        // Safe after validate_environment_templates above.
        let end = after_start
            .find('}')
            .expect("validated environment template");
        let name = &after_start[..end];
        let value = std::env::var(name)
            .map_err(|_| format!("required environment variable {name} is not set"))?;
        expanded.push_str(&value);
        remainder = &after_start[end + 1..];
    }
    expanded.push_str(remainder);
    Ok(expanded)
}
