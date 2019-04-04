use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::Strand;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::cmp::min;

#[derive(PartialEq)]
pub enum LibraryType {
    R1Sense,
    R2Sense,
    RandOrUnk,
}


pub trait LocusLike {
    fn as_contig(&self, use_strand: bool, as_frag: bool, sd: &ScaffoldDict, lt: LibraryType) -> Contig<String, Strand>;
}

impl LocusLike for bam::Record {
    fn as_contig(&self, use_strand: bool, frag: bool, sd: &ScaffoldDict, lt: LibraryType) -> Contig<String, Strand> {

        // only treat as frag is specified and possible
        let as_frag = match self.is_proper_pair() {
            true => frag,
            false => false,
        };

        // get the leftmost position of a template.
        let start = match as_frag {
            false => self.pos(),
            true  => min(self.pos(), self.mpos())
        };

        // get the rightmost pos of a frag or read,
        // either the cigar informed end of a single read
        // or the isize informed end of a pair.
        let end = match as_frag {
            false => self.cigar()
                         .end_pos()
                         .unwrap_or(self.seq().len() as i32 + start),
            true => start + self.insert_size().abs(),
        };

        let strand = match use_strand {
            true => match self.is_reverse() {
                true => Strand::Reverse,
                false => Strand::Forward,
            },
            false => Strand::Unknown,
        };

        let strand = if !use_strand {
            Strand::Unknown
        } else if !as_frag {
            match self.is_reverse() {
                true => Strand::Reverse,
                false => Strand::Forward
            }
        } else if lt == LibraryType::RandOrUnk {
            match self.is_reverse() {
                true => Strand::Reverse,
                false => Strand::Forward
            }
        } else if lt == LibraryType::R1Sense {
            if (self.is_first_in_template() & self.is_reverse()) {
                Strand::Reverse
            } else {
                Strand::Forward
            }
        } else if lt == LibraryType::R2Sense {
            if (self.is_first_in_template() & self.is_reverse()) {
                Strand::Forward
            } else {
                Strand::Reverse
            }
        } else {
            Strand::Unknown
        };


        Contig::new(
            sd.id_to_str(self.tid()).unwrap().to_owned(),
            start as isize,
            (end - start) as usize,
            strand,
        )
    }
}




#[cfg(test)]
mod tests {
    use crate::reads::locuslike::*;
    use crate::utility::scaffold_dict::ScaffoldDict;
    use bio::data_structures::annot_map::AnnotMap;
    use bio_types::annot::contig::Contig;
    use bio_types::strand::ReqStrand;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn nostrand_rec2contig_works() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let truth = "chr1:564476-564511";
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| a.as_contig(false, false, &tidmap, LibraryType::RandOrUnk));

        for i in res.into_iter() {
            assert_eq!(i.to_string(), truth);
        }
    }
}
