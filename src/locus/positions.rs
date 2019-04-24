use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::Pos;
use bio_types::strand;

/// Struct for iterator of positions within a Contig.
#[derive(Debug)]
pub struct Positions<'a> {
    region: &'a Contig<String,strand::ReqStrand>,
    idx: isize,
}

// This returns and Options<Result> for now, in case I want to implement more error checking later
/// Iterator for Pos within a Contig.
impl<'a> Iterator for Positions<'a> {
    type Item = Result<Pos<String, strand::ReqStrand>, &'static str>;
    fn next(&mut self) -> Option<Result<Pos<String, strand::ReqStrand>, &'static str>>
        {
            let w = self.region.length();
        match (self.idx as usize) < w {
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
    use crate::locus::from_rec::*;
    use crate::utility::scaffold_dict::ScaffoldDict;
    use bio_types::annot::contig::Contig;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    use crate::locus::positions::*;


    #[test]
    fn positions_from_contig() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(2)
            .map(|a| Contig::from_read(&a, false, &tidmap))
            .map(|a| a.unwrap());

        for i in res.into_iter() {
            println!("\nREAD: {:?}",i);
            i.positions().for_each(|a| println!("{:?}",a));
        }
    }
}
