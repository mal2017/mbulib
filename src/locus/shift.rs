use bio_types::annot::pos::Pos;
use bio_types::annot::loc::*;
use bio_types::annot::contig::Contig;
use bio_types::strand::*;

/// Shift a Loc implementing object. Always relative to strand, unless
/// Loc object is unstranded in which case it is treated as Forward
/// stranded by convention.
pub trait Shift : Loc {
    fn shift(&self, d: isize) -> Self;
}


impl Shift for Contig<String, ReqStrand> {
    fn shift(&self, d: isize) -> Self {
        let d2 = match self.strand() {
            ReqStrand::Reverse => -1 * d,
            _ => d,
        };
        Contig::new(self.refid().to_string(), self.start() + d2, self.length(), self.strand())
    }
}


impl Shift for Pos<String, ReqStrand> {
    fn shift(&self, d: isize) -> Self {
        let d2 = match self.strand() {
            ReqStrand::Reverse => -1 * d,
            _ => d,
        };
        Pos::new(self.refid().to_string(), self.start() + d2, self.strand())
    }
}



#[cfg(test)]
mod tests {
    use crate::locus::from_rec::*;
    use crate::locus::shift::*;
    use crate::utility::scaffold_dict::ScaffoldDict;
    use bio::data_structures::annot_map::AnnotMap;
    use bio_types::annot::contig::Contig;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn nostrand_shift_contig() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let truth = "chr1:564481-564516(+)";
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| Contig::from_read(&a, false, &tidmap))
            .map(|a| a.shift(5));

        for i in res.into_iter() {
            assert_eq!(i.to_string(), truth);
        }

    }
}
