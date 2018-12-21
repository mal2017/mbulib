use rust_htslib::bam::{Read, Records};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::ReadError;
use std::iter::Iterator;

#[derive(Debug)]
pub struct Sort<R: Iterator> {
    records: R,
}


impl<R: Iterator> Iterator for Sort<R> {
    type Item = <R as Iterator>::Item;

    #[inline]
    fn next(&mut self) -> Option<<R as Iterator>::Item> { self.records.next() }
}



pub trait Sortable : Iterator + Sized {
    fn by_read_name(self) -> Sort<Self>;
}


impl<'a, R: 'a + Read> Sortable for Records<'a, R>  {
    fn by_read_name(self) -> Sort<Self> {
        // TODO make this a sorting iterator
        Sort {
            records: self
        }
    }
}







#[cfg(test)]
mod tests {
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    use crate::bam::sort::*;

    #[test]
    fn it_works() {
        let bampath = Path::new("test/aln.mini.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let z = bam.records()
                   .by_read_name()
                   .map(|a| {println!("{:?}",a); a})
                   .count();
        assert_eq!(2 + 2, 4);
    }
}
