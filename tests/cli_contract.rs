use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_percolator-rs"))
}

fn fixture() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/sample.pin")
}

fn run(arguments: &[&str]) -> Output {
    binary().args(arguments).output().unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).unwrap()
}

#[test]
fn no_input_reports_the_input_contract_and_exit_code() {
    let output = run(&[]);
    assert_eq!(output.status.code(), Some(2));
    let message = stderr(&output);
    assert!(message.starts_with("usage: percolator-rs [flags] input.pin"));
    assert!(message.contains("Separate target/decoy searches (mix-max) are not supported."));
}

#[test]
fn explicit_options_override_profiles_in_either_order() {
    for arguments in [
        vec!["--maxiter", "1", "--fast", fixture()],
        vec!["--fast", "--maxiter", "1", fixture()],
    ] {
        let output = run(&arguments);
        assert!(output.status.success(), "{}", stderr(&output));
        assert!(stderr(&output)
            .contains("profile: fast (model=svm, maxiter=1, subset-max-train=20000)"));
    }
}

#[test]
fn legacy_svm_model_aliases_remain_compatible() {
    for alias in ["svm", "linear"] {
        let output = run(&["--model", alias, "--maxiter", "1", fixture()]);
        assert!(output.status.success(), "{}", stderr(&output));
        assert!(stderr(&output).contains("profile: canonical (model=svm, maxiter=1"));
    }
}

#[test]
fn removed_neural_options_fail_clearly() {
    for (arguments, expected) in [
        (
            vec!["--model", "mlp", fixture()],
            "unknown --rescore-model 'mlp' (only svm is available)",
        ),
        (
            vec!["--model", "neural", fixture()],
            "unknown --rescore-model 'neural' (only svm is available)",
        ),
        (
            vec!["--mlp-hidden", "8", fixture()],
            "--mlp-hidden is no longer supported; neural rescoring has been removed",
        ),
    ] {
        let output = run(&arguments);
        assert_eq!(output.status.code(), Some(2));
        assert_eq!(stderr(&output).trim(), expected);
    }
}

#[test]
fn invalid_scientific_option_values_keep_their_diagnostics() {
    let cases = [
        (
            vec!["--protein-inference", "guess", fixture()],
            "unknown --protein-inference 'guess' (use picked|bayesian)",
        ),
        (
            vec!["--null-target-win-prob", "0", fixture()],
            "invalid --null-target-win-prob (must be finite and in (0, 1))",
        ),
        (
            vec!["--model", "guess", fixture()],
            "unknown --rescore-model 'guess' (only svm is available)",
        ),
    ];
    for (arguments, expected) in cases {
        let output = run(&arguments);
        assert_eq!(output.status.code(), Some(2));
        assert_eq!(stderr(&output).trim(), expected);
    }
}

#[test]
fn mutually_exclusive_modes_fail_before_input_loading() {
    let output = run(&["--join", "--ensemble", "missing-a.pin", "missing-b.pin"]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        stderr(&output).trim(),
        "--ensemble and --join are mutually exclusive"
    );
}

#[test]
fn malformed_ensemble_input_preserves_its_contract() {
    let output = run(&["--ensemble", "invalid", "engine=also-missing.pin"]);
    assert_eq!(output.status.code(), Some(2));
    let message = stderr(&output);
    assert!(message.contains("error: invalid ensemble input 'invalid'; use ENGINE=PIN"));
    assert!(!message.contains("parse error"));
}

#[test]
fn help_and_version_are_successful_without_an_input() {
    for option in ["--help", "-h"] {
        let output = run(&[option]);
        assert!(output.status.success());
        let help = String::from_utf8(output.stdout).unwrap();
        assert!(help.contains("Usage: percolator-rs"));
        assert!(help.contains("--results-psms"));
        assert!(help.contains("--protein-inference"));
        assert!(output.stderr.is_empty());
    }
    let output = run(&["--version"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap().trim(),
        concat!("percolator-rs ", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn invalid_options_cannot_silently_change_the_analysis() {
    for arguments in [
        vec!["--typo", fixture()],
        vec!["--seed", "oops", fixture()],
        vec!["--maxiter", "oops", fixture()],
        vec!["--maxiter", "-1", fixture()],
        vec!["--subset-max-train", "oops", fixture()],
        vec!["--num-threads", "oops", fixture()],
        vec!["--num-threads", "0", fixture()],
        vec!["--cpos", "oops", fixture()],
        vec!["--cpos", "-1", fixture()],
        vec!["--cneg", "0", fixture()],
        vec!["--cpos", "NaN", fixture()],
        vec!["--cneg", "inf", fixture()],
        vec!["--results-psms"],
        vec!["--results-psms", "--fast", fixture()],
    ] {
        let output = run(&arguments);
        assert_eq!(output.status.code(), Some(2), "{arguments:?}");
        assert!(!stderr(&output).contains("parsed "), "{arguments:?}");
        assert!(!stderr(&output).contains("panicked"), "{arguments:?}");
    }
}

#[test]
fn output_paths_cannot_overwrite_inputs_or_other_results() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.pin");
    let original = std::fs::read(fixture()).unwrap();
    std::fs::write(&input, &original).unwrap();
    let same = input.to_str().unwrap();
    let output = run(&["--results-psms", same, same]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("refers to an input file"));
    assert_eq!(std::fs::read(&input).unwrap(), original);

    let result = directory.path().join("results.tsv");
    let result = result.to_str().unwrap();
    let output = run(&["--results-psms", result, "--results-peptides", result, same]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("more than one output"));
    assert!(!std::path::Path::new(result).exists());

    let alias = directory.path().join("alias.pin");
    std::fs::hard_link(&input, &alias).unwrap();
    let output = run(&["--results-psms", alias.to_str().unwrap(), same]);
    #[cfg(unix)]
    assert_eq!(output.status.code(), Some(2));
    #[cfg(unix)]
    assert_eq!(std::fs::read(&input).unwrap(), original);
}

#[test]
fn invalid_output_directory_has_an_actionable_error() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("missing").join("results.tsv");
    let output = run(&["--results-psms", path.to_str().unwrap(), fixture()]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("invalid output path"));
    assert!(!stderr(&output).contains("panicked"));
}

#[cfg(target_os = "linux")]
#[test]
fn write_failure_is_reported_without_aborting() {
    for option in [
        "--results-psms",
        "--results-peptides",
        "--results-proteins",
        "--feature-report",
    ] {
        let output = run(&["--maxiter", "0", option, "/dev/full", fixture()]);
        assert_eq!(output.status.code(), Some(1), "{option}");
        assert!(stderr(&output).contains("error (/dev/full)"));
        assert!(!stderr(&output).contains("panicked"));
    }
}

#[test]
fn joined_inputs_require_matching_feature_names_and_order() {
    let directory = tempfile::tempdir().unwrap();
    let a = directory.path().join("a.pin");
    let b = directory.path().join("b.pin");
    for (path, features) in [(&a, "score\tother"), (&b, "other\tscore")] {
        std::fs::write(path, format!(
            "SpecId\tLabel\tScanNr\t{features}\tPeptide\tProteins\nrow\t1\t1\t2\t3\tK.PEPTIDE.R\tP1\n"
        )).unwrap();
    }
    let output = run(&["--join", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("identical feature names in the same order"));
    assert!(!stderr(&output).contains("panicked"));
}

#[test]
fn joined_per_file_counts_describe_the_reported_psms() {
    let directory = tempfile::tempdir().unwrap();
    let first = directory.path().join("first.pin");
    let second = directory.path().join("second.pin");
    let results = directory.path().join("targets.tsv");
    std::fs::copy(fixture(), &first).unwrap();
    std::fs::copy(fixture(), &second).unwrap();
    let output = run(&[
        "--join",
        "--maxiter",
        "1",
        "--results-psms",
        results.to_str().unwrap(),
        first.to_str().unwrap(),
        second.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{}", stderr(&output));
    let per_source: usize = stderr(&output)
        .lines()
        .filter(|line| line.starts_with("  ["))
        .map(|line| line.rsplit_once(' ').unwrap().1.parse::<usize>().unwrap())
        .sum();
    let text = std::fs::read_to_string(results).unwrap();
    let accepted = text
        .lines()
        .skip(1)
        .filter(|line| line.split('\t').nth(2).unwrap().parse::<f64>().unwrap() < 0.01)
        .count();
    assert!(accepted > 0);
    assert_eq!(per_source, accepted);
}
