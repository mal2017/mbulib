use bio::io::bed;
use rust_htslib::bam;
use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::annot::loc::*;
use bio_types::annot::pos::Pos;
use bio_types::strand;
use rust_htslib::bam::Records;
use rust_htslib::bam::HeaderView;
use crate::locus::from_read::*;
use crate::locus::shift::*;
use crate::utility::scaffold_dict::ScaffoldDict;
use rust_htslib::bam::Read;
use std::cmp::{min, max};
use crate::locus::shift::Shift;

// #[derive(Debug)]
// pub struct Positions<'a, R: 'a + Loc>
//     where
//         R::RefID: Clone,
//         R::Strand: Copy + Into<strand::Strand>,
//     {
//     region: &'a R,
//     idx: isize,
//     pos: &'a Pos<String, strand::Strand>,
// }
//
// impl<'a, L: Loc> Iterator for Positions<'a, L>
//     where
//         L::RefID: Clone,
//         L::Strand: Copy + Into<strand::Strand>,
//     {
//     type Item = Result<Pos<String, strand::Strand>, &'static str>;
//
//     fn next(&mut self) -> Option<Result<Pos<String, strand::Strand>, &'static str>>
//         {
//         let sp = self.region.start();
//
//         let len = self.length()s;
//
//         //let p = self.region.first_pos();
//         self.idx += 1;
//         Some(Ok(pos))
//
//     }
// }

#[derive(Debug)]
pub struct Positions<'a> {
    region: &'a Contig<String,strand::ReqStrand>,
    idx: isize,
}

impl<'a> Iterator for Positions<'a> {
    type Item = Result<Pos<String, strand::Strand>, &'static str>;
    fn next(&mut self) -> Option<Result<Pos<String, strand::Strand>, &'static str>>
        {
        let sp = self.region.start();
        let w = self.region.length();

        let strand: strand::Strand = self.region.strand().into();

        let refid: String = self.region.refid().to_string();

        let pos = Pos::new(refid, self.idx, strand);
        self.idx += 1;
        Some(Ok(pos))

    }
}

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
    use crate::locus::from_read::*;
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
        let truth = "chr1:564476-564511(+)";
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| Contig::from_read(&a, false, &tidmap));

        for i in res.into_iter() {
            i.positions().for_each(|a| println!("{:?}",a));
        }
    }
}
