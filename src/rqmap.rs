use rust_htslib::bam;
use bio::data_structures::annot_map::*;
use bio_types::annot::contig::Contig;
use bio_types::annot::pos::Pos;
use bio_types::annot::loc::Loc;
use bio_types::strand::ReqStrand;
use rust_htslib::bam::HeaderView;
use crate::locus::from_rec::*;
use crate::scaffold_dict::ScaffoldDict;
use crate::locus::positions::PositionScan;
use rust_htslib::bam::Read;
use std::cmp::{min, max};
use crate::errors::*;

trait AppendRecord {
    fn append(&mut self, r: &bam::Record, as_frags: bool, sd: &ScaffoldDict, pf: &Option<fn(Contig<String,ReqStrand>)-> Contig<String,ReqStrand>>) -> Result<(), InvalidRecordError>;
}


impl AppendRecord for AnnotMap<String, Contig<String,ReqStrand>> {
    fn append(&mut self, r: &bam::Record, as_frags: bool, sd: &ScaffoldDict, pf: &Option<fn(Contig<String,ReqStrand>)-> Contig<String,ReqStrand>>) -> Result<(), InvalidRecordError> {

        match pf {
            None => self.insert_loc(Contig::from_read(r, as_frags, sd)?),
            Some(f) => self.insert_loc(f(Contig::from_read(r, as_frags, sd)?)),
        }

        Ok(())
    }
}

#[derive(Debug)]
pub enum LibraryType {
    R1Sense,
    R2Sense,
    Unstranded,
    ATAC,
    DNASE,

}

#[derive(Debug)]
pub enum CountStrand {
    Plus,
    Minus,
    Both,
}

/// Struct holds a library of NGS reads as an AnnotMap
/// and strandedness information.
#[derive(Debug)]
pub struct RQMap {
    construction: LibraryType,
    plus: AnnotMap<String, Contig<String,ReqStrand>>,
    minus: AnnotMap<String, Contig<String,ReqStrand>>
}

// TODO: clean up this with enum logic
impl RQMap {
    /// Retrieve the counts within a locus for reads mapped to both strands.
    pub fn counts_within<D: Into<ReqStrand>>(&self, p: &Contig<String, D>, ct: &CountStrand) -> usize {
        match ct {
            CountStrand::Plus => self.plus.find(p).count(),
            CountStrand::Minus => self.minus.find(p).count(),
            CountStrand::Both => self.plus.find(p).count() + self.minus.find(p).count(),
        }
    }

    /// Retrieve the coverage at a pointlike locus for reads mapped to both strands.
    pub fn coverage_at<D: Into<ReqStrand>>(&self, p: &Pos<String, D>, ct: &CountStrand) -> usize {
        match ct {
            CountStrand::Plus => self.plus.find(p).count(),
            CountStrand::Minus => self.minus.find(p).count(),
            CountStrand::Both => self.plus.find(p).count() + self.minus.find(p).count(),
        }
    }

    pub fn quantify<L>(&self, p: &L, ct: &CountStrand) -> usize
        where
            L: Loc<RefID = String>,
    {
        match ct {
            CountStrand::Plus => self.plus.find(p).count(),
            CountStrand::Minus => self.minus.find(p).count(),
            CountStrand::Both => self.plus.find(p).count() + self.minus.find(p).count(),
        }
    }

    /// Retrieve the coverage across a contiguous locus for reads mapped to both strands.
    pub fn coverage_across(&self, c: &Contig<String, ReqStrand>, ct: &CountStrand) -> Vec<usize> {
        c.positions()
         .map(|a| a.unwrap())
         .map(|a| self.coverage_at(&a, ct))
         .collect()
    }

    /// Create an RQMap from a bam reader.
    // TODO: descriptive error on unmapped?
    pub fn from_reader(mut b: bam::Reader,
                          as_frags: bool,
                          lt: LibraryType,
                          rf: Option<fn(&bam::Record) -> bool>,
                          pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Result<Self,InvalidRecordError> {
        let mut plus: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let mut minus: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let hd: HeaderView = b.header().clone();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(&hd);
        let mut r: bam::Record = bam::Record::new();

        // map
        while let Ok(_r) = b.read(&mut r) {
            match rf {
                None => { match r.is_reverse() {
                    true => minus.append(&r, as_frags, &sd, &pf)?,
                    false => plus.append(&r, as_frags, &sd, &pf)?,
                } },
                Some(f) => {
                    match f(&r) {
                        true => { match r.is_reverse() {
                            true => minus.append(&r, as_frags, &sd, &pf)?,
                            false => plus.append(&r, as_frags, &sd, &pf)?,
                        } },
                        false => continue,
                    }
                }
            };

        }

        Ok(RQMap {
            construction: lt,
            plus: plus,
            minus: minus
        })
    }

    /// Create an RQMap from an indexed bam reader.
    // TODO: descriptive error on unmapped?
    pub fn from_indexed(mut b: bam::IndexedReader,
                            as_frags: bool,
                            c: Vec<Contig<String,ReqStrand>>,
                            lt: LibraryType,
                            rf: Option<fn(&bam::Record) -> bool>,
                            pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Result<Self,InvalidRecordError> {
        let mut plus: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let mut minus: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let hd: HeaderView = b.header().clone();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(&hd);

        let mut chr: u32;
        let mut c1: u32;
        let mut c2: u32;
        let mut r: bam::Record = bam::Record::new();

        for x in c.into_iter() {
            chr = match sd.str_to_id(&x.refid()) {
                Some(i) => i as u32,
                None => continue,
            };
            c1 = x.first_pos().start() as u32;
            c2 = x.last_pos().start() as u32;

            match b.fetch(chr, min(c1,c2), max(c1,c2)) {
                Ok(_x) => {} ,
                _ => return Err(InvalidRecordError::NoneSuchRecord),
            };

            while let Ok(_r) = b.read(&mut r) {
                match rf {
                    None => { match r.is_reverse() {
                        true => minus.append(&r, as_frags, &sd, &pf)?,
                        false => plus.append(&r, as_frags, &sd, &pf)?,
                    } },
                    Some(f) => {
                        match f(&r) {
                            true => { match r.is_reverse() {
                                true => minus.append(&r, as_frags, &sd, &pf)?,
                                false => plus.append(&r, as_frags, &sd, &pf)?,
                            } },
                            false => continue,
                        }
                    }
                };

            }
        }

        Ok(RQMap {
            construction: lt,
            plus: plus,
            minus: minus,
        })
    }

}


#[cfg(test)]
mod tests {
    use bio::data_structures::annot_map::AnnotMap;
    use rust_htslib::bam;
    use std::path::Path;
    use rust_htslib::bam::Read;
    use bio_types::annot::contig::Contig;
    use bio_types::annot::loc::Loc;
    use crate::rqmap::*;
    use crate::locus::shift::*;
    use rust_htslib::bam::HeaderView;
    use bio_types::strand::ReqStrand;

    fn tn5shift(c: Contig<String,ReqStrand>) -> Contig<String,ReqStrand> {
        let new = match c.strand() {
            ReqStrand::Forward => c.shift(4),
            ReqStrand::Reverse => c.shift(5),
        };
        new.first_pos().contig()
    }

    fn mapq_filt(b: &bam::Record) -> bool {
        b.mapq() > 30
    }

    #[test]
    fn rqmap_from_reader() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::Reader::from_path(bampath).unwrap();

        // TODO Work on this test
        let _r = RQMap::from_reader(bam, false, LibraryType::Unstranded, None, None);

    }

    #[test]
    fn rqmap_from_reader_filt() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::Reader::from_path(bampath).unwrap();


        // TODO work on this test
        let _r = RQMap::from_reader(bam, false, LibraryType::Unstranded, Some(mapq_filt), None);

    }
    #[test]
    fn rqmap_from_reader_preproc() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::Reader::from_path(bampath).unwrap();


        // TODO work on this test
        let _r = RQMap::from_reader(bam, false, LibraryType::Unstranded, None, Some(tn5shift));

    }

    #[test]
    fn rqmap_from_indexed() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::IndexedReader::from_path(bampath).unwrap();


        let c1: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                     1000000,
                                                     1000000,
                                                     ReqStrand::Forward);

        let _r = RQMap::from_indexed(bam,
                                         false,
                                         vec!(c1),
                                         LibraryType::Unstranded,
                                         None,
                                         None);
        }

    #[test]
    fn rqmap_from_indexed_filt() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::IndexedReader::from_path(bampath).unwrap();


        let c1: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                         1000000,
                                                         1000000,
                                                         ReqStrand::Forward);

        let _r = RQMap::from_indexed(bam,
                                             false,
                                             vec!(c1),
                                             LibraryType::Unstranded,
                                             Some(mapq_filt),
                                             None);
    }

    #[test]
    fn rqmap_from_indexed_preproc() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::IndexedReader::from_path(bampath).unwrap();


        let c1: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                         1000000,
                                                         1000000,
                                                         ReqStrand::Forward);

        let _r = RQMap::from_indexed(bam,
                                             false,
                                             vec!(c1),
                                             LibraryType::Unstranded,
                                             None,
                                             Some(tn5shift));
    }

    #[test]
    fn coverage_across_region() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::IndexedReader::from_path(bampath).unwrap();

        let c0: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                     564475,
                                                     60,
                                                     ReqStrand::Forward);

        let r = RQMap::from_indexed(bam,
                                        false,
                                        vec!(c0.clone()),
                                        LibraryType::Unstranded,
                                        None,
                                        None);

        }

}
