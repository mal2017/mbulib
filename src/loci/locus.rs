use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::Strand;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::bed;


/// LocusLike wraps various record types that might be reasonably
/// represented as a range. For example: rust_htslib::bam::Record.
/// This allows us to generate generics for many different
// TODO make this a reference??
pub struct LocusLike<T> {
    r: T
}

/// Construct a LocusLike struct from a bam record.
impl LocusLike<bam::Record> {
    pub fn from_bam_rec(r: bam::Record) -> Self {
        Self {
            r: r,
        }
    }
    pub fn as_contig(&self, use_strand: bool, sd: &ScaffoldDict) -> Contig<String, Strand> {
        let start = self.r.pos();
        let end = self
            .r
            .cigar()
            .end_pos()
            .unwrap_or(self.r.seq().len() as i32 + start);
        let w = end - start;

        let strand = match use_strand {
            true => match self.r.is_reverse() {
                true => Strand::Reverse,
                false => Strand::Forward,
            },
            false => Strand::Unknown,
        };

        Contig::new(
            sd.id_to_str(self.r.tid()).unwrap().to_owned(),
            start as isize,
            w as usize,
            strand,
        )
    }
}

/// Construct a LocusLike struct from a bed record.
impl LocusLike<bed::Record> {
    pub fn from_bed_rec(r: bed::Record) -> Self {
        Self {
            r: r,
        }
    }
    pub fn as_contig(&self, use_strand: bool) -> Contig<String, Strand> {

        let strand = match use_strand {
            true => match self.r.strand() {
                Some(s) => s,
                None => Strand::Unknown,
            },
            false => Strand::Unknown,
        };

        Contig::new(
            self.r.chrom().to_string(),
            self.r.start() as isize,
            (self.r.end() - self.r.start()) as usize,
            strand
        )
    }

}


#[cfg(test)]
mod tests {
    use crate::loci::locus::*;
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
            .map(|a| LocusLike::from_bam_rec(a))
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
            .map(|a| LocusLike::from_bam_rec(a))
            .map(|a| a.as_contig(false, &tidmap)).collect();
    }
}
