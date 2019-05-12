// Copyright 2018-2019 Matt Lawlor.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.
//!
//! RQMap is a crate for reading bam files into a queryable map object (Read Query Map).
//! Methods are provided for generating RQMaps, filtering or modifying reads, and quantifying
//! signal within regions.

extern crate bio_types;
extern crate rust_htslib;

#[macro_use]
extern crate quick_error;
//extern crate lazysort

/// Utilities for common housekeeping tasks, such as
/// file header reading/maniplation.
pub mod scaffold_dict;

/// Utilities for manipulating structs that implement 'loc'.
pub mod locus;

/// Utilities for converting bams to usable data structures.
pub mod rqmap;

/// Useful and descriptive errors
pub mod errors;
