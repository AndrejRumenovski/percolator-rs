//! Import externally prepared PINs with an explicit reproducible recipe.
//! Recipes describe work already performed; they are never shell-executed.
use super::{
    cache::{Cache, Object},
    download::{hash_file, Budgets},
    *,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::Path,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Recipe {
    pub steps: Vec<Lineage>,
}
#[derive(Debug, Serialize)]
pub struct PreparationPlan {
    pub accession: Pxd,
    pub input_bytes: u64,
    pub additional_cache_bytes: u64,
    pub download_bytes: u64,
    pub cache_current_bytes: u64,
    pub cache_limit_bytes: u64,
    pub free_filesystem_bytes: u64,
    pub peak_working_bytes: u64,
    pub expected_evictions: Vec<String>,
    pub retained: String,
}
pub fn validate_recipe(m: &Manifest, recipe: &Recipe) -> Result<()> {
    if recipe.steps.is_empty() {
        return Err("preparation recipe must contain at least one step".into());
    }
    let mut known: BTreeSet<_> = m
        .inventory
        .iter()
        .map(|f| f.id.clone())
        .chain(m.lineage.iter().map(|s| s.id.clone()))
        .collect();
    for s in &recipe.steps {
        if s.id.is_empty()
            || known.contains(&s.id)
            || s.inputs.is_empty()
            || s.inputs.iter().any(|i| !known.contains(i))
        {
            return Err("recipe step IDs must be unique and reference remote file IDs or preceding lineage steps".into());
        }
        if s.tool.is_empty()
            || s.tool_version.as_ref().is_none_or(String::is_empty)
            || s.parameters.is_empty()
        {
            return Err(
                "each preparation step needs an explicit tool, version and reproduction parameters"
                    .into(),
            );
        }
        if s.kind == "database_search"
            && (s.protein_database.is_none()
                || s.database_sha256
                    .as_ref()
                    .is_none_or(|s| s.len() != 64 || !s.bytes().all(|b| b.is_ascii_hexdigit()))
                || s.decoy_generation.is_none())
        {
            return Err("database_search recipe requires protein database, SHA-256 and decoy generation method".into());
        }
        known.insert(s.id.clone());
    }
    if recipe.steps.last().unwrap().kind != "pin" {
        return Err("last preparation step must have kind 'pin'".into());
    }
    Ok(())
}
pub fn import_plan(
    cache: &Cache,
    m: &Manifest,
    pin: &Path,
    recipe: &Recipe,
    b: &Budgets,
) -> Result<PreparationPlan> {
    validate_recipe(m, recipe)?;
    if fs::canonicalize(pin)?.starts_with(&cache.root) {
        return Err("import expects an external PIN; managed objects are already tracked".into());
    }
    if !fs::metadata(pin)?.is_file() {
        return Err("prepared PIN must be a regular file".into());
    }
    let size = fs::metadata(pin)?.len();
    let peak = total([size, b.safety])?;
    if size > cache.limit || b.max_working_space.is_some_and(|n| peak > n) {
        return Err("prepared PIN import exceeds cache/working-space budget".into());
    }
    let s = cache.status()?;
    let evictions = cache.eviction_plan(size, peak, &BTreeSet::new())?;
    Ok(PreparationPlan{accession:m.accession.clone(),input_bytes:size,additional_cache_bytes:size,download_bytes:0,cache_current_bytes:s.large_data_bytes,cache_limit_bytes:cache.limit,free_filesystem_bytes:s.free_filesystem_bytes,peak_working_bytes:peak,expected_evictions:evictions,retained:"Imported PIN follows retention; recipe, hashes and remote ancestry survive eviction. The user-supplied input file is not modified.".into()})
}
pub fn import(
    cache: &mut Cache,
    m: &mut Manifest,
    pin: &Path,
    recipe: Recipe,
    b: &Budgets,
    retention: Retention,
) -> Result<String> {
    cache.require_write()?;
    import_plan(cache, m, pin, &recipe, b)?;
    m.preparation_attempts.push(PreparationAttempt {
        steps: recipe.steps.clone(),
        state: "running".into(),
        error: None,
        started_unix_seconds: now(),
    });
    cache.save_manifest(m)?;
    let result = import_inner(cache, m, pin, recipe, b, retention);
    let attempt = m.preparation_attempts.last_mut().unwrap();
    match &result {
        Ok(_) => attempt.state = "verified".into(),
        Err(error) => {
            attempt.state = "failed".into();
            attempt.error = Some(error.to_string());
        }
    }
    cache.save_manifest(m)?;
    result
}
fn import_inner(
    cache: &mut Cache,
    m: &mut Manifest,
    pin: &Path,
    recipe: Recipe,
    b: &Budgets,
    retention: Retention,
) -> Result<String> {
    cache.require_write()?;
    let plan = import_plan(cache, m, pin, &recipe, b)?;
    drop(crate::pin::parse(pin.to_str().ok_or("non-UTF8 PIN path")?)?);
    let before = hash_file(pin)?;
    let key = format!("prepared-{}", before.sha256);
    // Protect a possible existing identical prepared artifact while making room.
    let evictions = cache.eviction_plan(
        before.bytes,
        total([before.bytes, b.safety])?,
        &BTreeSet::from([key.clone()]),
    )?;
    cache.evict(&evictions, false)?;
    let id = recipe.steps.last().unwrap().id.clone();
    let final_rel = format!("prepared/{key}");
    let dest = cache.path(&final_rel)?;
    let mut o = Object {
        key: key.clone(),
        relative_path: format!("tmp/{key}.part"),
        bytes: before.bytes,
        state: State::Partial,
        local_sha256: None,
        verification: vec![],
        projects: BTreeSet::from([m.accession.clone()]),
        last_used_unix_seconds: now(),
        retention,
        result_verified: false,
        reproducible: true,
        etag: None,
        pin_validated: false,
    };
    if let Some(old) = cache.index.objects.get(&key) {
        o.projects.extend(old.projects.clone());
        if old.retention == Retention::Keep {
            o.retention = Retention::Keep;
        }
    }
    if dest.exists() {
        if hash_file(&dest)?.sha256 != before.sha256 {
            return Err("existing prepared object is corrupt".into());
        }
    } else {
        cache.index.objects.insert(key.clone(), o.clone());
        cache.save_index()?;
        let part = cache.path(&o.relative_path)?;
        let mut input = File::open(pin)?;
        let mut output = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&part)?;
        let mut buf = [0u8; 65536];
        let mut written = 0u64;
        loop {
            let n = input.read(&mut buf)?;
            if n == 0 {
                break;
            }
            if total([written, n as u64])? > plan.input_bytes {
                return Err("prepared input changed/grew during import".into());
            }
            if cache::available_space(&cache.root)? < total([n as u64, b.safety])? {
                return Err("free-space safety margin reached during PIN import".into());
            }
            output.write_all(&buf[..n])?;
            written += n as u64;
        }
        output.sync_all()?;
        if hash_file(&part)?.sha256 != before.sha256 {
            return Err("prepared input changed during import; partial retained".into());
        }
        fs::rename(part, &dest)?;
        File::open(dest.parent().unwrap())?.sync_all()?;
    }
    o.relative_path = final_rel;
    o.local_sha256 = Some(before.sha256.clone());
    o.state = State::Prepared;
    o.pin_validated = true;
    m.local_files.insert(id.clone(), o.local());
    cache.index.objects.insert(key.clone(), o);
    cache.save_index()?;
    let mut steps = recipe.steps;
    steps.last_mut().unwrap().output_sha256 = Some(before.sha256.clone());
    m.lineage.extend(steps);
    m.prepared_pins.insert(
        id.clone(),
        PreparedPin {
            id: id.clone(),
            object_key: key,
            sha256: before.sha256,
            bytes: before.bytes,
            lineage_id: id.clone(),
            retention,
        },
    );
    cache.save_manifest(m)?;
    Ok(id)
}
pub fn run_plan(
    cache: &Cache,
    m: &Manifest,
    id: &str,
    b: &Budgets,
    result_bytes: u64,
) -> Result<PreparationPlan> {
    let pin = m
        .prepared_pins
        .get(id)
        .ok_or("unknown prepared PIN ID; inspect the manifest")?;
    let o = cache
        .index
        .objects
        .get(&pin.object_key)
        .ok_or("prepared PIN absent; regenerate from retained recipe")?;
    if !matches!(o.state, State::Prepared) || !cache.path(&o.relative_path)?.is_file() {
        return Err(
            "prepared PIN evicted/corrupt; regenerate and import with a new recipe ID".into(),
        );
    }
    let peak = total([pin.bytes, result_bytes, b.safety])?;
    if b.max_working_space.is_some_and(|n| peak > n) {
        return Err("prepared analysis exceeds working-space budget".into());
    }
    let s = cache.status()?;
    let evictions = cache.eviction_plan(
        result_bytes,
        total([result_bytes, b.safety])?,
        &BTreeSet::from([pin.object_key.clone()]),
    )?;
    Ok(PreparationPlan {
        accession: m.accession.clone(),
        input_bytes: pin.bytes,
        additional_cache_bytes: 0,
        download_bytes: 0,
        cache_current_bytes: s.large_data_bytes,
        cache_limit_bytes: cache.limit,
        free_filesystem_bytes: s.free_filesystem_bytes,
        peak_working_bytes: peak,
        expected_evictions: evictions,
        retained: "Verified final results and recipe; PIN follows its retention and pinning policy"
            .into(),
    })
}
