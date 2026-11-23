// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Error types for vep-io.

/// Errors that can occur during I/O operations.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IoError {
    /// An underlying I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// A record a parser cannot read, or a value a writer cannot serialise, the
    /// reason as the whole message: it names the field refused, and a refused
    /// line is skipped with it under every input format.
    #[error("{0}")]
    VcfParse(String),

    /// Input refused by design, with the reason as the whole message.
    #[error("{0}")]
    Unsupported(String),

    /// A transcript cache an input names a record in but a lookup cannot read,
    /// the message the parser's text followed by the cache's own refusal; the
    /// run stops and prints it under the `cache` component, since the input
    /// itself is sound.
    #[error("{0}")]
    CacheUnreadable(String),

    /// A line Ensembl VEP's parser warns about and reads past, the warning's
    /// text as the whole message.
    #[error("{0}")]
    SkippedLine(String),

    /// A Parquet file that cannot be written or read, with the encoder's or
    /// decoder's reason.
    #[error("Parquet: {0}")]
    Parquet(String),
}
