extern crate bio_types;
extern crate itertools;
extern crate linear_map;
extern crate rust_htslib;
//extern crate lazysort

/// Utilities for common housekeeping tasks, such as
/// file header reading/maniplation.
pub mod utility;

/// Utilities for manipulating structs that implement 'loc'.
pub mod locus;

/// Utilities for converting bams to usable data structures.
pub mod ngslibrary;

/// Useful and descriptive errors
pub mod errors;
