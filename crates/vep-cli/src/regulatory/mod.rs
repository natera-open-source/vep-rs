// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Regulatory overlaps (Ensembl VEP 116.3): the regulatory tree of a converted
//! cache and its loader (`source`), the GFF files `--regulatory_gff` reads in
//! its place (`gff`), the per-chromosome overlap index (`index`) and the rows
//! a variant gains (`rows`).

mod gff;
mod index;
mod rows;
mod source;

pub use gff::CELL_TYPE_WARNING;
pub use rows::CellTypes;
pub use source::{
    count_shard_features, BindingMatrices, Feature, FeatureKind, RegulatorySource,
    GFF_REPLACES_CACHE_WARNING,
};

/// The regulatory fixture cache the module's tests share.
#[cfg(test)]
pub(crate) mod fixtures {
    use std::path::Path;
    use std::sync::Arc;

    use vep_core::cache_info::REGULATORY_SCHEMA;
    use vep_core::cache_manifest::CacheManifest;
    use vep_core::variant::InputVariant;

    use super::RegulatorySource;
    use crate::cache_resolve::ResolvedCache;

    /// The fixtures' cell-type list, the names the fixture shard's activity
    /// string is aligned to.
    pub fn cell_types() -> Vec<String> {
        serde_json::from_slice(include_bytes!(
            "../../tests/fixtures/regulatory/cell_types.json"
        ))
        .unwrap()
    }

    /// A cache at `dir` whose metadata declares regulatory data in the current
    /// shape, with the fixtures' cell-type list.
    pub fn resolved_cache(dir: &Path) -> ResolvedCache {
        ResolvedCache {
            dir: dir.to_string_lossy().into_owned(),
            info: {
                let mut info = vep_core::cache_info::CacheInfo::default();
                info.regulatory = true;
                info.regulatory_schema = Some(REGULATORY_SCHEMA);
                info.cell_types = cell_types();
                info
            },
            manifest: Arc::new(CacheManifest::load(dir).unwrap()),
            assembly: Some("GRCh38".into()),
            fasta: None,
            synonyms: Arc::new(Default::default()),
            source_type: crate::cache_resolve::SourceType::Ensembl,
            notices: Vec::new(),
        }
    }

    /// Writes the regulatory tree under `dir`: the fixtures' matrices table,
    /// the hand-built shard on chromosome 21 as `1-1000000.json`, and `extra`
    /// shards (name, body) beside it.
    pub fn write_fixture_cache(dir: &Path, extra: &[(&str, &str)]) {
        let chr = dir.join("regulatory/21");
        std::fs::create_dir_all(&chr).unwrap();
        std::fs::write(
            dir.join("regulatory/binding_matrices.json"),
            include_bytes!("../../tests/fixtures/regulatory/binding_matrices.json"),
        )
        .unwrap();
        std::fs::write(
            chr.join("1-1000000.json"),
            include_bytes!("../../tests/fixtures/regulatory/1-1000000.json"),
        )
        .unwrap();
        for (name, body) in extra {
            std::fs::write(chr.join(name), body).unwrap();
        }
    }

    /// The fixture cache with its shard written twice, as two overlapping
    /// shards, so the loader must emit each feature once.
    pub fn source() -> (tempfile::TempDir, Arc<RegulatorySource>) {
        let dir = tempfile::tempdir().unwrap();
        let fixture = include_str!("../../tests/fixtures/regulatory/1-1000000.json");
        write_fixture_cache(dir.path(), &[("1000001-2000000.json", fixture)]);
        let source = RegulatorySource::open(&resolved_cache(dir.path()), true)
            .unwrap()
            .unwrap();
        (dir, source)
    }

    pub fn variant(start: u64, end: u64, reference: &[u8], alternate: &[u8]) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            reference.to_vec(),
            alternate.to_vec(),
        )
    }
}
