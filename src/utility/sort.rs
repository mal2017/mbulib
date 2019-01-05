use itertools::Itertools;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{Read, Records};
use std::cmp::Ordering;
use std::iter::Iterator;
use std::vec::IntoIter;
//use mbulib::RecordSort;

/// Trait for sorting Record Iterators various ways.
pub trait RecordSort {
    type Item;

    /// Sort Records by name, ID, etc. See implementations for
    /// specific structs.
    fn name_sort(self) -> IntoIter<Self::Item>;
}

/// Implementations for sort iterators over bam Records.
impl<'a, R: 'a + Read> RecordSort for Records<'a, R> {
    type Item = Record;

    /// Sort a Records object by read name.
    fn name_sort(self) -> IntoIter<Record> {
        self.map(|a| a.unwrap())
            .sorted_by(|a, b| cmp_read_name(&a, &b))
            .into_iter()
    }
}

// Internal method for comparing Bam Record vs Bam Record lexicographically
// by read name.
fn cmp_read_name(rec1: &Record, rec2: &Record) -> Ordering {
    let read1name = String::from_utf8(rec1.qname().to_owned()).unwrap();
    let read2name = String::from_utf8(rec2.qname().to_owned()).unwrap();
    read1name.cmp(&read2name)
}

#[cfg(test)]
mod tests {
    use crate::utility::sort::*;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn name_sort_bam_works() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let res = bam
            .records()
            .into_iter()
            .name_sort()
            .map(|a| String::from_utf8(a.qname().to_owned()).unwrap());

        for (x, y) in res.tuple_windows() {
            assert!(y >= x);
        }
    }
}
