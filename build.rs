// Build identity is provenance only; it never changes scientific configuration.
fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    if let Ok(output) = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
    {
        if output.status.success() {
            let commit = String::from_utf8_lossy(&output.stdout);
            let commit = commit.trim();
            if commit.len() == 40 && commit.bytes().all(|c| c.is_ascii_hexdigit()) {
                println!("cargo:rustc-env=PERCOLATOR_RS_BUILD_COMMIT={commit}");
            }
        }
    }
}
