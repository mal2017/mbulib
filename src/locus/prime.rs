use bio_types::annot::pos::Pos;
use bio_types::annot::loc::*;
use bio_types::strand::*;

/*
// define EitherEnd for all types implementing loc
pub trait EitherEnd : Loc {
    /// Find the 5` end of a region.
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy;

    /// Find the 3` end of a region.
    fn threeprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy;
}


impl<T: Loc> EitherEnd for T {
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy,
    {
        match self.strand().into() {
            ReqStrand::Reverse => self.last_pos(),
            _ => self.first_pos(),

        }
    }

    fn threeprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy,
    {
        match self.strand().into() {
            ReqStrand::Reverse => self.first_pos(),
            _ => self.last_pos(),

        }
    }
}
*/

// define EitherEnd for all types implementing loc
pub trait EitherEnd : Loc {
    /// Find the 5` end of a region.
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy;

}


impl<T: Loc> EitherEnd for T {
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Copy,
        Self::Strand: Into<ReqStrand> + Copy,
    {
        match self.strand().into() {
            ReqStrand::Reverse => self.last_pos(),
            _ => self.first_pos(),

        }
    }

}

#[cfg(test)]
mod tests {
    use crate::reads::locuslike::*;
    use crate::locus::prime::*;
    use crate::utility::scaffold_dict::ScaffoldDict;
    use bio::data_structures::annot_map::AnnotMap;
    use bio_types::annot::contig::Contig;
    use bio_types::strand::ReqStrand;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn nostrand_shift_contig() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let truth = "chr1:564476";
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| a.as_contig(false, false, &tidmap, LibraryType::RandOrUnk))
            .map(|a| a.fiveprime());

        for i in res.into_iter() {
            assert_eq!(i.to_string(), truth);
        }

    }
}
