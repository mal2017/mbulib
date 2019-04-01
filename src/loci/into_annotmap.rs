use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::HeaderView;
use crate::loci::locus::LocusLike;
use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::Strand;

// for reference:
// https://stackoverflow.com/questions/30630810/using-generic-iterators-instead-of-specific-list-types
// https://github.com/mal2017/mbulib/commit/690a67d4be92dc3fb3dff8588d955a0e3c7723bb

pub trait IntoAnnotMap {
    fn collect_annotmap(self) -> AnnotMap<String, String> ;
}

// impl<T> IntoAnnotMap for T
//     where T: Iterator<Item = bam::Record>
//     {
//     fn collect_annotmap(self) -> AnnotMap<String, String> {
//         let mut map: AnnotMap<String,String> = AnnotMap::new();
//         map
//     }
// }

impl IntoAnnotMap for rust_htslib::bam::Reader {
    fn collect_annotmap(mut self) -> AnnotMap<String, String> {
        let mut map: AnnotMap<String,String> = AnnotMap::new();
        let hd: &HeaderView = self.header();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);

        let z = self.records()
                    .map(|a| a.unwrap())
                    .map(|a| LocusLike::from_bam_rec(a))
                    .map(|a| a.as_contig(true, &sd));


        map
    }
}




use bio::io::bed;

#[cfg(test)]
mod tests {
    use crate::loci::locus::*;
    use crate::loci::into_annotmap::IntoAnnotMap;
    use crate::utility::scaffold_dict::ScaffoldDict;
    use bio::data_structures::annot_map::AnnotMap;
    use bio_types::annot::contig::Contig;
    use bio_types::strand::ReqStrand;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    #[test]
    fn bam_reads_2_annotmap() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();

        let res = bam
            .collect_annotmap();


    }
}
