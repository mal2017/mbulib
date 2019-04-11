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
pub struct Positions {
    region: Contig<String,strand::Strand>,
    idx: isize,
}

impl Iterator for Positions {
    type Item = Result<Pos<String, strand::Strand>, &'static str>;
    fn next(&mut self) -> Option<Result<Pos<String, strand::Strand>, &'static str>>
        {
        let sp = self.region.start();

        let strand: strand::Strand = self.region.strand().into();

        

        let pos = Pos::new("chr1".to_string(), 100, strand::Strand::Unknown);
        self.idx += 1;
        Some(Ok(pos))

    }
}
