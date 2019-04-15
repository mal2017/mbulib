use bio::io::bed;
use rust_htslib::bam;
use bio::data_structures::annot_map::*;
use bio_types::annot::contig::Contig;
use bio_types::annot::pos::Pos;
use bio_types::annot::loc::Loc;
use bio_types::strand::{ReqStrand, Strand};
use rust_htslib::bam::Records;
use rust_htslib::bam::HeaderView;
use crate::locus::from_read::*;
use crate::utility::scaffold_dict::ScaffoldDict;
use crate::locus::positions::PositionScan;
use rust_htslib::bam::Read;
use std::cmp::{min, max};

trait IntoAnnotMap{
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>, hd: &HeaderView, pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Result<(), &'static str>;
}

/// TODO: ADD OPTIONAL FILTERING/PREPROC STEP FOR RAW READS
impl<T: Iterator<Item=bam::Record>> IntoAnnotMap for T {
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>, hd: &HeaderView, pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Result<(), &'static str> {
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);

        // TODO make this return meaningful error
        match pf {
            None => {self.map(|a| Contig::from_read(&a, false, &sd)).for_each(|a| am.insert_loc(a));},
            Some(f) => {self.map(|a| Contig::from_read(&a, false, &sd)).map(|a| f(a)).for_each(|a| am.insert_loc(a));},
        }

        Ok(())
    }
}

#[derive(Debug)]
pub enum LibraryType {
    R1Sense,
    R2Sense,
    Unstranded,
}

/// Struct holds a library of NGS reads as an AnnotMap, as well as the original reader struct
/// and strandedness information.
#[derive(Debug)]
pub struct NGSLibrary<T>
    where T: bam::Read {
    reader: T,
    construction: LibraryType,
    map: AnnotMap<String, Contig<String,ReqStrand>>,
}

impl<T: Read> NGSLibrary<T> {
    /// Retrieve the coverage at a pointlike locus.
    pub fn coverage_at<D: Into<ReqStrand>>(&self, p: &Pos<String, D>) -> usize {
        self.map.find(p).count()
    }
    /// Retrieve the coverage across a contiguous locus.
    pub fn coverage_across(&self, c: &Contig<String, ReqStrand>) -> Vec<usize> {
        c.positions()
         .map(|a| a.unwrap())
         .map(|a| self.coverage_at(&a))
         .collect()
    }
}


impl NGSLibrary<bam::Reader> {
    /// Create an NGSLibrary from an indexed bam reader.
    pub fn from_reader(mut b: bam::Reader, lt: LibraryType, pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Self {
        let mut map: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let hd: HeaderView = b.header().clone();
        let _r = b.records()
                    .map(|a| a.unwrap())
                    .add_to(&mut map, &hd, pf);
        NGSLibrary {
            reader: b,
            construction: lt,
            map: map,
        }
    }
}

impl NGSLibrary<bam::IndexedReader> {
    /// Create an NGSLibrary from an indexed bam reader.
    pub fn from_indexed(mut b: bam::IndexedReader, c: Vec<Contig<String,ReqStrand>>, lt: LibraryType, pf: Option<fn(Contig<String,ReqStrand>) -> Contig<String,ReqStrand>>) -> Self {
        let mut map: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let hd: HeaderView = b.header().clone();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(&hd);

        for x in c.into_iter() {
            let chr = sd.str_to_id(&x.refid()).unwrap() as u32;
            let c1 = x.first_pos().start() as u32;
            let c2 = x.last_pos().start() as u32;
            b.fetch(chr, min(c1,c2), max(c1,c2));
            b.records()
             .map(|a| a.unwrap())
             .add_to(&mut map, &hd, pf);
        }

        NGSLibrary {
            reader: b,
            construction: lt,
            map: map,
        }
    }
}




#[cfg(test)]
mod tests {
    use bio::data_structures::annot_map::AnnotMap;
    use bio::io::bed;
    use rust_htslib::bam;
    use std::path::Path;
    use rust_htslib::bam::Read;
    use bio_types::annot::contig::Contig;
    use bio_types::annot::loc::Loc;
    use bio_types::annot::contig::*;
    use crate::ngslibrary::*;
    use crate::locus::shift::*;
    use rust_htslib::bam::HeaderView;
    use bio_types::strand::ReqStrand;

    #[test]
    fn reads_into_annotmap() {
        let mut map: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hd: HeaderView = bam.header().clone();
        let res = bam
            .records()
            .take(5)
            .map(|a| a.unwrap())
            .add_to(&mut map, &hd, None);

        assert_eq!(res, Ok(()))
    }

    #[test]
    fn reads_into_ngslib() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();

        let r = NGSLibrary::from_reader(bam, LibraryType::Unstranded, None);

    }

    #[test]
    fn reads_into_ngslib_preproc() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();

        fn tn5shift(c: Contig<String,ReqStrand>) -> Contig<String,ReqStrand> {
            match c.strand() {
                ReqStrand::Forward => c.shift(4),
                ReqStrand::Reverse => c.shift(-5),
            }

        }

        let r = NGSLibrary::from_reader(bam, LibraryType::Unstranded, Some(tn5shift));

    }

    #[test]
    fn indexed_reads_into_ngslib() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::IndexedReader::from_path(bampath).unwrap();

        let c0: Contig<String,ReqStrand> = Contig::new("chr10".to_string(),
                                                     1000000,
                                                     1000000,
                                                     ReqStrand::Forward);

        let c1: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                     1000000,
                                                     1000000,
                                                     ReqStrand::Forward);

        let r = NGSLibrary::from_indexed(bam,
                                         vec!(c0,c1),
                                         LibraryType::Unstranded,
                                         None);
        }

    #[test]
    fn coverage_across_region() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::IndexedReader::from_path(bampath).unwrap();

        let c0: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                     564442,
                                                     5000,
                                                     ReqStrand::Forward);

        let r = NGSLibrary::from_indexed(bam,
                                        vec!(c0.clone()),
                                        LibraryType::Unstranded,
                                        None);
        println!("COV ACROSS: {:?}", r);

        println!("COV ACROSS: {:?}", r.coverage_across(&c0));

        }
}
