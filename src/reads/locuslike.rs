use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::Strand;
use rust_htslib::bam;
use rust_htslib::bam::Read;

pub trait LocusLike {
    fn as_contig(&self, use_strand: bool, sd: &ScaffoldDict) -> Contig<String, Strand>;
}

impl LocusLike for bam::Record {
    fn as_contig(&self, use_strand: bool, sd: &ScaffoldDict) -> Contig<String, Strand> {
        let start = self.pos();
        let end = self
            .cigar()
            .end_pos()
            .unwrap_or(self.seq().len() as i32 + start);
        let w = end - start;

        let strand = match use_strand {
            true => match self.is_reverse() {
                true => Strand::Reverse,
                false => Strand::Forward,
            },
            false => Strand::Unknown,
        };

        Contig::new(
            sd.id_to_str(self.tid()).unwrap().to_owned(),
            start as isize,
            w as usize,
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
    fn bam_rec_2_contig_works() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let truth = "chr1:564476-564511(+)";
        let res = bam
            .records()
            .map(|a| a.unwrap())
            .take(1)
            .map(|a| a.as_contig(true, &tidmap));

        for i in res.into_iter() {
            assert_eq!(i.to_string(), truth);
        }
    }

    #[test]
    fn make_locus_works() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let hdrv = bam.header();
        let tidmap = ScaffoldDict::from_header_view(&hdrv);
        let res: Vec<_> = bam
            .records()
            .map(|a| a.unwrap())
            .take(2)
            .map(|a| a.as_contig(false, &tidmap)).collect();
    }
}
