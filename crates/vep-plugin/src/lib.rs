// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin trait definition, dylib loader, and registry for VEP.
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

/// Feature types that a plugin can annotate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureType {
    Transcript,
    RegulatoryFeature,
    MotifFeature,
}

impl FeatureType {
    /// The feature types a plugin declares, from the comma-separated list its
    /// `vep_plugin_feature_types` symbol returns: the words Ensembl VEP's
    /// `feature_types` takes (`BaseVepPlugin.pm` `_check_types`), where
    /// `Feature` is every feature class and `Intergenic` names the rows that
    /// have no feature, which the host hands to no plugin. Any other word is a
    /// load error.
    pub fn parse_list(text: &str) -> Result<Vec<Self>, PluginError> {
        let mut types = Vec::new();
        for word in text.split(',').map(str::trim).filter(|w| !w.is_empty()) {
            let declared: &[Self] = match word {
                "Transcript" => &[Self::Transcript],
                "RegulatoryFeature" => &[Self::RegulatoryFeature],
                "MotifFeature" => &[Self::MotifFeature],
                "Feature" => &[
                    Self::Transcript,
                    Self::RegulatoryFeature,
                    Self::MotifFeature,
                ],
                "Intergenic" => &[],
                other => {
                    return Err(PluginError::Load(format!(
                        "unknown feature type {other:?} in vep_plugin_feature_types"
                    )))
                }
            };
            for kind in declared {
                if !types.contains(kind) {
                    types.push(*kind);
                }
            }
        }
        Ok(types)
    }

    /// The row kind this feature type annotates.
    pub fn row_kind(self) -> vep_core::consequence::FeatureType {
        match self {
            Self::Transcript => vep_core::consequence::FeatureType::Transcript,
            Self::RegulatoryFeature => vep_core::consequence::FeatureType::RegulatoryFeature,
            Self::MotifFeature => vep_core::consequence::FeatureType::MotifFeature,
        }
    }
}

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
