// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Error types for vep-io.

/// Errors that can occur during I/O operations.
#[derive(Debug, thiserror::Error)]
pub enum IoError {
    /// An underlying I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// A VCF parsing error.
    #[error("VCF parse error: {0}")]
    VcfParse(String),

    /// Input refused by design, with the reason as the whole message.
    #[error("{0}")]
    Unsupported(String),

    /// A line Ensembl VEP's parser warns about and reads past, the warning's
    /// text as the whole message.
    #[error("{0}")]
    SkippedLine(String),
}
