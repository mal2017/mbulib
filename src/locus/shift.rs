use bio_types::annot::pos::Pos;
use bio_types::annot::loc::*;
use bio_types::annot::contig::Contig;
use bio_types::strand::*;

// Shift a Loc implementing object. Always relative to strand, unless
// Loc object is unstranded in which case it is treated as Forward
// stranded by convention.
pub trait Shift : Loc {
    fn shift(self, d: isize) -> Self;
}



impl Shift for Contig<String, Strand> {
    fn shift(self, d: isize) -> Self {
        let d2 = match self.strand() {
            Strand::Reverse => -1 * d,
            _ => d,
        };
        Contig::new(self.refid().to_string(), self.start() + d2, self.length(), self.strand())
    }
}


impl Shift for Pos<String, Strand> {
    fn shift(self, d: isize) -> Self {
        let d2 = match self.strand() {
            Strand::Reverse => -1 * d,
            _ => d,
        };
        Pos::new(self.refid().to_string(), self.start() + d2, self.strand())
    }
}
