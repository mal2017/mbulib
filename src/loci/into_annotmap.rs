use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::HeaderView;
use crate::loci::locus::LocusLike;
use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::strand::Strand;

// for reference:
// https://stackoverflow.com/questions/30630810/using-generic-iterators-instead-of-specific-list-types
// https://github.com/mal2017/mbulib/commit/690a67d4be92dc3fb3dff8588d955a0e3c7723bb

pub trait IntoAnnotMap {
    fn collect_annotmap(self) -> AnnotMap<String, Contig<String, Strand>> ;
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,Strand>>) -> Result<(), &'static str>;
}

impl<T: Read> IntoAnnotMap for T {
    fn collect_annotmap(mut self) -> AnnotMap<String, Contig<String,Strand>> {
        let mut map: AnnotMap<String,Contig<String,Strand>> = AnnotMap::new();
        let hd: &HeaderView = self.header();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);

        let z = self.records()
                    .map(|a| a.unwrap())
                    .map(|a| a.as_contig(true, &sd));

        for i in z.into_iter() {
            map.insert_loc(i);
        }

        map
    }
    fn add_to(mut self, am: &mut AnnotMap<String, Contig<String,Strand>>) -> Result<(), &'static str> {
        let hd: &HeaderView = self.header();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);
        let z = self.records()
                    .map(|a| a.unwrap())
                    .map(|a| a.as_contig(true, &sd));

        for i in z.into_iter() {
            am.insert_loc(i);
        }
        Ok(())
    }
}


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
        //let bampath = Path::new("test/hs.pe.test.bam");
        let bampath = Path::new("../test.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();

        let res = bam
            .collect_annotmap();


    }
}
