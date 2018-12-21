use rust_htslib::bam::{Read, Records};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::ReadError;

#[derive(Debug)]
struct Sort<'a, R: 'a + Iterator> {
    records: &'a mut R,
}

impl<'a, R: Iterator> Iterator for Sort<'a, R> {
    //type Item = Result<Record, ReadError>;
    type Item = <R as Iterator>::Item;

    #[inline]
    fn next(&mut self) -> Option<<R as Iterator>::Item> { self.records.next() }
}


/*
pub trait Sortable {
    fn by_read_name(self) -> Sort<Self, F>;
}


impl<'a, R: Read> Sortable for Records<'a, R> {
    fn sort_by_name<B, F>(self) -> Sort<Self, F> where {

        self
    }
}
*/
