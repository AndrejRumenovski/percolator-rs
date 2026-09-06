//! Public PRIDE research infrastructure. No scientific algorithms live here.
pub mod cache;
pub mod client;
pub mod download;
pub mod model;
pub mod prepare;
pub mod workflow;
pub use model::*;
pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

pub fn now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

pub fn bytes(value: &str) -> std::result::Result<u64, String> {
    let value = value.trim().to_ascii_uppercase();
    let end = value
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(value.len());
    let n: u64 = value[..end]
        .parse()
        .map_err(|_| "use an integer byte count, e.g. 20GB or 512MiB")?;
    let factor = match &value[end..] {
        "" | "B" => 1,
        "KB" => 1_000,
        "MB" => 1_000_000,
        "GB" => 1_000_000_000,
        "TB" => 1_000_000_000_000,
        "KIB" => 1 << 10,
        "MIB" => 1 << 20,
        "GIB" => 1 << 30,
        "TIB" => 1 << 40,
        _ => return Err("unknown size unit; use B, MB, GB, TB, MiB or GiB".into()),
    };
    n.checked_mul(factor)
        .ok_or_else(|| "size overflows u64".into())
}

pub fn total(values: impl IntoIterator<Item = u64>) -> Result<u64> {
    values.into_iter().try_fold(0u64, |sum, v| {
        sum.checked_add(v)
            .ok_or_else(|| "storage size overflow".into())
    })
}
