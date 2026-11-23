// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The features of a genome, read from its converted Ensembl VEP cache: the genetic code
//! of every coding transcript, the biotypes, the coding sequences cut off at an end, the
//! sequences, the synonyms and the variant, SIFT and regulatory data the cache carries.
//! [`extract`] reads one cache into a [`features::Features`], [`shard`] reads what the
//! cache reader of `vep-seeded-variants` leaves out of a transcript, and [`gather`]
//! writes one table row a genome from a directory of such files.

pub mod extract;
pub mod features;
pub mod gather;
pub mod shard;
