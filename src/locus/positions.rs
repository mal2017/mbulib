use bio_types::strand;
use bio_types::annot::{contig::Contig, loc::Loc, pos::Pos};

/// This struct holds state for an iterator over the positions within a contiguous region.
#[derive(Debug, Clone)]
pub struct Positions<'a> {
    region: &'a Contig<String,strand::ReqStrand>,
    idx: isize,
}

/// This iterator produces the positions within a contiguous region.
impl<'a> Iterator for Positions<'a> {
    type Item = Result<Pos<String, strand::ReqStrand>, &'static str>;
    fn next(&mut self) -> Option<Result<Pos<String, strand::ReqStrand>, &'static str>>
        {
        match (self.idx as usize) < self.region.length() {
            true => {
                let strand: strand::ReqStrand = self.region.strand().into();
                let refid: String = self.region.refid().to_string();
                let pos_in = Pos::new(refid, self.idx, strand);
                let pos = self.region.pos_outof(&pos_in);
                self.idx += 1;
                Some(Ok(pos.unwrap())) // TODO unscrew up this please, make it match and return error/option intelligently
            },
            false => {
                None
            }
        }
    }
}

/// Trait for generating an iterator of positions.
pub trait PositionScan {
    fn positions(&self) -> Positions;
}

impl PositionScan for Contig<String, strand::ReqStrand> {
    fn positions(&self) -> Positions {
        Positions {
            region: self,
            idx: 0
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{locus::from_rec::LocFromRec, scaffold_dict, locus::positions::*};
    use bio_types::annot::contig::Contig;
    use rust_htslib::bam::*;
    use std::path::Path;
    use crate::locus::positions::*;

    #[test]
    fn positions_from_contig() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = scaffold_dict::ScaffoldDict::from_header_view(&hdrv);
        let res: Vec< Contig<String,strand::ReqStrand> > = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| Contig::from_read(&a, false, &tidmap))
            .map(|a| a.unwrap()).collect();



    }
}
