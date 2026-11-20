// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The C ABI a VEP plugin dylib exports, the dylib loader and the registry of
//! loaded plugins, with [`VepPlugin`], the in-process plugin trait.
//!
//! Plugins are compiled as cdylib crates and loaded at runtime.
//! The [`abi`] module defines the C FFI contract, [`loader`] handles dylib
//! loading, and [`registry`] manages the set of active plugins.

pub mod abi;
pub mod loader;
pub mod registry;

pub use loader::LoadedPlugin;
pub use registry::PluginRegistry;

use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// The kind of row a plugin runs on: vep-core's `RowKind`, which the built-in
/// plugins' `feature_types` return as well.
pub use vep_core::consequence::RowKind as FeatureType;

/// Errors from plugin execution.
#[derive(Debug, thiserror::Error)]
pub enum PluginError {
    #[error("plugin initialization failed: {0}")]
    Init(String),
    #[error("plugin execution failed: {0}")]
    Run(String),
    #[error("plugin load failed: {0}")]
    Load(String),
    #[error("incompatible plugin ABI version: expected {expected}, got {got}")]
    AbiMismatch { expected: u32, got: u32 },
}

/// In-process plugin interface; loaded plugins export the C ABI in [`abi`], built-ins use `BuiltinPlugin`.
pub trait VepPlugin: Send + Sync {
    fn name(&self) -> &str;
    fn version(&self) -> &str;
    fn feature_types(&self) -> &[FeatureType];
    fn header_info(&self) -> Vec<(String, String)>;
    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<Vec<(String, String)>, PluginError>;
    fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
        Ok(())
    }
    fn finish(&mut self) -> Result<(), PluginError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_type_names_round_trip_and_reject_any_other_string() {
        for kind in [
            FeatureType::Transcript,
            FeatureType::RegulatoryFeature,
            FeatureType::MotifFeature,
            FeatureType::Intergenic,
        ] {
            assert_eq!(FeatureType::parse(kind.as_str()), Some(kind));
        }
        assert_eq!(FeatureType::parse("transcript"), None);
        assert_eq!(FeatureType::parse("Gene"), None);
    }

    #[test]
    fn a_consequence_feature_type_maps_to_the_plugin_kind_of_the_same_name() {
        use vep_core::consequence::FeatureType as Core;
        for core in [
            Core::Transcript,
            Core::RegulatoryFeature,
            Core::MotifFeature,
        ] {
            assert_eq!(FeatureType::from(core).as_str(), core.as_str());
        }
    }
}
