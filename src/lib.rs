extern crate bio_types;
extern crate itertools;
extern crate linear_map;
extern crate rust_htslib;
//extern crate lazysort

/// Utilities for common housekeeping tasks, such as
/// file header reading/maniplation.
pub mod utility;

/// High-level utilities for importing/exporting genomic data types.
pub mod io;

/// Mid-level utilities for dealing with converting records Into
/// locus-like data structures
pub mod loci;
