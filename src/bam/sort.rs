use rust_htslib::bam::{Read, Records};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::ReadError;
use std::iter::Iterator;
use std::cmp::Ordering;
use std::vec::IntoIter;
use itertools::Itertools;


pub trait BamSort {
    type Item;

    fn name_sort(self) -> IntoIter<Self::Item>;
}

impl<'a, R: 'a + Read> BamSort for Records<'a, R> {
    type Item = Record;

    fn name_sort(self) -> IntoIter<Record> {
        self.map(|a| a.unwrap())
            .sorted_by(|a ,b| {cmp_read_name(&a, &b)} )
    }
}


fn cmp_read_name(rec1: &Record, rec2: &Record) -> Ordering {
    let read1name = String::from_utf8(rec1.qname().to_owned()).unwrap();
    let read2name = String::from_utf8(rec2.qname().to_owned()).unwrap();
    read1name.cmp(&read2name)
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    use crate::bam::sort::*;
    use rust_htslib::prelude::*;
    use std::str;
    #[test]
    fn it_works() {
        let bampath = Path::new("test/aln.mini.bam");
        let opath = Path::new("/Users/mlawlor/Desktop/test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdr = bam::header::Header::from_template(bam.header());
        let hdrv = bam.header().to_owned();
        let mut obam = bam::Writer::from_path(opath, &hdr).unwrap();
        bam.records().into_iter()
                     .name_sort()
                     .map(|a| println!("{:?}", String::from_utf8(a.qname().to_owned()).unwrap())).for_each(drop);
                     //.map(|a| obam.write(&a).unwrap()).for_each(drop);
        assert_eq!(2 + 2, 4);
    }
}
