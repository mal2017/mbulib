//use bio_types::annot::contig::Contig;
use bio_types::annot::pos::Pos;
use bio_types::annot::loc::*;
use bio_types::strand::*;

// define EitherEnd for all types implementing loc
pub trait EitherEnd : Loc {
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>;
    //fn threeprime(self) -> Self;
}

impl<T: Loc> EitherEnd for T {
    fn fiveprime(&self) -> Pos<Self::RefID, Self::Strand>
    //where
    //    Self::RefID: Clone,
    //    Self::Strand: Into<ReqStrand> + Copy
    {
        match self.strand().into() {
            ReqStrand::Forward => self.first_pos(),
            ReqStrand::Reverse => self.last_pos(),
        }

    }
}
