use percolator_rs::benchmark_manifest::{DatasetRegistry, ManifestError};

const VALID: &str = r#"
version = 1

[[datasets]]
id = "PXD032157"
pride_accession = "PXD032157"
source = "PRIDE Archive"
organism = "Anopheles gambiae metaproteome"
experiment_type = "DDA"
instrument = "Q Exactive HF"
search_engine = "Comet"
pin_path = "${PERCOLATOR_BENCH_DATA}/PXD032157/**/*.pin"
file_count = 65
approximate_input_size = "2.30 GiB"
protein_level_evaluation = false
notes = "Protein evaluation is unsuitable for this metaproteomic search."
preparation = "Place data outside the repository."
"#;

#[test]
fn parses_the_committed_pxd032157_registry() {
    let registry =
        DatasetRegistry::load("bench/datasets.toml").expect("committed manifest is valid");
    assert_eq!(registry.version, 1);
    assert_eq!(registry.datasets.len(), 1);
    assert_eq!(registry.datasets[0].id, "PXD032157");
    assert_eq!(registry.datasets[0].file_count, Some(65));
    assert!(!registry.datasets[0].protein_level_evaluation);
}

#[test]
fn parses_all_supported_metadata() {
    let registry = DatasetRegistry::from_toml(VALID).expect("valid manifest");
    assert_eq!(
        registry.datasets[0].pride_accession.as_deref(),
        Some("PXD032157")
    );
    assert_eq!(
        registry.datasets[0].instrument.as_deref(),
        Some("Q Exactive HF")
    );
}

#[test]
fn reports_multiple_validation_errors() {
    let invalid = VALID
        .replace("id = \"PXD032157\"", "id = \"bad id\"")
        .replace("file_count = 65", "file_count = 0")
        .replace("${PERCOLATOR_BENCH_DATA}", "${bad-name}");
    let error = DatasetRegistry::from_toml(&invalid).expect_err("manifest should fail validation");
    let message = error.to_string();
    assert!(matches!(error, ManifestError::Validation(_)));
    assert!(message.contains("id must use ASCII"));
    assert!(message.contains("file_count must be greater than zero"));
    assert!(message.contains("invalid ${...} environment template"));
}

#[test]
fn rejects_missing_required_fields_and_unknown_keys() {
    let missing = VALID.replace("source = \"PRIDE Archive\"\n", "");
    let error = DatasetRegistry::from_toml(&missing).expect_err("source is required");
    assert!(error.to_string().contains("missing field `source`"));

    let unknown = format!("{VALID}\nnot_a_dataset_field = true\n");
    let error = DatasetRegistry::from_toml(&unknown).expect_err("unknown fields should fail");
    assert!(error
        .to_string()
        .contains("unknown field `not_a_dataset_field`"));
}
