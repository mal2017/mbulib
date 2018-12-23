use itertools::Itertools;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::ReadError;
use rust_htslib::bam::{Read, Records};
use std::cmp::Ordering;
use std::iter::Iterator;
use std::vec::IntoIter;

/// Trait for sorting Records various ways.
pub trait BamSort {
    type Item;

    fn name_sort(self) -> IntoIter<Self::Item>;
}

impl<'a, R: 'a + Read> BamSort for Records<'a, R> {
    type Item = Record;

    /// Sort a Records object by read name.
    fn name_sort(self) -> IntoIter<Record> {
        self.map(|a| a.unwrap())
            .sorted_by(|a, b| cmp_read_name(&a, &b))
    }
}

// Internal method for comparing Record vs Record lexicographically
// by read name.
fn cmp_read_name(rec1: &Record, rec2: &Record) -> Ordering {
    let read1name = String::from_utf8(rec1.qname().to_owned()).unwrap();
    let read2name = String::from_utf8(rec2.qname().to_owned()).unwrap();
    read1name.cmp(&read2name)
}

#[cfg(test)]
mod tests {
    use crate::bam::sort::*;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn name_sort_bam_works() {
        let bampath = Path::new("test/aln.mini.bam");
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
