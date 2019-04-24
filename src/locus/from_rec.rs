use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::strand::ReqStrand;
use rust_htslib::bam;
use std::cmp::min;


pub trait LocFromRec {
    /// From a bam record create a struct implementing @Loc.
    fn from_read(rec: &bam::Record, frag: bool, sd: &ScaffoldDict) -> Option<Self>
    where
        Self: Loc + Sized;
}

impl LocFromRec for Contig<String, ReqStrand> {
    fn from_read(rec: &bam::Record, frag: bool, sd: &ScaffoldDict) -> Option<Self>  {
        if rec.is_unmapped() {
            return None;
        }

        // only treat as frag is specified and possible
        let as_frag = match rec.is_proper_pair() {
            true => frag,
            false => false,
        };

        // get the leftmost position of a template.
        let start = match as_frag {
            false => rec.pos(),
            true  => min(rec.pos(), rec.mpos())
        };

        // get the rightmost pos of a frag or read,
        // either the cigar informed end of a single read
        // or the isize informed end of a pair.
        let end = match as_frag {
            false => rec.cigar()
                         .end_pos()
                         .unwrap_or(rec.seq().len() as i32 + start),
            true => start + rec.insert_size().abs(),
        };

        let strand = match rec.is_reverse() {
                true => ReqStrand::Reverse,
                false => ReqStrand::Forward,
        };

        Some(Contig::new(
            sd.id_to_str(rec.tid()).unwrap().to_owned(),
            start as isize,
            (end - start) as usize,
            strand,
        ))
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum InvalidRecordError {
        UnmappedRecord {
            description("Record is unmapped/unpaired")
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


    #[test]
    fn contig_from_read() {
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
            assert_eq!(i.unwrap().to_string(), truth);
        }
    }
}
