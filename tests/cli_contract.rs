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
        assert!(stderr(&output).contains(
            "profile: fast (model=svm, maxiter=1, subset-max-train=20000)"
        ));
    }
}

#[test]
fn documented_model_aliases_select_the_same_model() {
    for alias in ["svm", "linear"] {
        let output = run(&["--model", alias, "--maxiter", "1", fixture()]);
        assert!(output.status.success(), "{}", stderr(&output));
        assert!(stderr(&output).contains("profile: canonical (model=svm, maxiter=1"));
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
            "unknown --rescore-model 'guess' (use svm|mlp)",
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
