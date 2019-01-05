extern crate itertools;
extern crate linear_map;
extern crate rust_htslib;
//extern crate lazysort


use std::vec::IntoIter;

/// Trait for sorting various Record Iterators various ways.
pub trait RecordSort {
    type Item;

    /// Sort Records by name, ID, etc.
    fn name_sort(self) -> IntoIter<Self::Item>;
}




pub mod bam;
