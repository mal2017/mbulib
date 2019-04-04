use bio_types::annot::pos::Pos;
use bio_types::annot::loc::*;
use bio_types::strand::*;

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
