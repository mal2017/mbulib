use rust_htslib::bam;
use bio::data_structures::annot_map::*;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::strand::ReqStrand;
use rust_htslib::bam::HeaderView;
use crate::locus::from_rec::*;
use crate::scaffold_dict::ScaffoldDict;
use crate::locus::positions::PositionScan;
use rust_htslib::bam::Read;
use std::cmp::{min, max};
use crate::errors::*;
use crate::library_strategy::*;
use crate::rqmap::RQMap;

pub trait AppendRecord {
    fn append(&mut self,
        r: &bam::Record,
        as_frags: bool,
        sd: &ScaffoldDict,
        pf: &Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) ->
        Result<(), InvalidRecordError>;
}


impl AppendRecord for AnnotMap<String, Contig<String,ReqStrand>> {
    fn append(&mut self,
        r: &bam::Record,
        as_frags: bool,
        sd: &ScaffoldDict,
        pf: &Option<fn(Contig<String,ReqStrand>)-> Contig<String,ReqStrand>>) ->
        Result<(), InvalidRecordError> {

        match pf {
            None => self.insert_loc(Contig::from_read(r, as_frags, sd)?),
            Some(f) => self.insert_loc(f(Contig::from_read(r, as_frags, sd)?)),
        }

        Ok(())
    }
}

impl AppendRecord for RQMap {
    fn append(&mut self,
        r: &bam::Record,
        as_frags: bool,
        sd: &ScaffoldDict,
        pf: &Option<fn(Contig<String,ReqStrand>)-> Contig<String,ReqStrand>>) ->
        Result<(), InvalidRecordError> {

        match pf {
            None => self.insert_loc(Contig::from_read(r, as_frags, sd)?),
            Some(f) => self.insert_loc(f(Contig::from_read(r, as_frags, sd)?)),
        }

        Ok(())
    }
}
