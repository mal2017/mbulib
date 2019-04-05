use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use rust_htslib::bam::Records;
use rust_htslib::bam::HeaderView;
use crate::locus::from_read::*;
use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;
use rust_htslib::bam::Read;

pub trait IntoAnnotMap{
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>, hd: &HeaderView) -> Result<(), &'static str>;
}

impl<T: Iterator<Item=bam::Record>> IntoAnnotMap for T {
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>, hd: &HeaderView) -> Result<(), &'static str> {
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);
        let z = self.map(|a| Contig::from_read(&a, false, &sd));

        for i in z.into_iter() {
            // for some reason this method doesn't return
            // a result.
            am.insert_loc(i);
        }
        Ok(())
    }
}



#[cfg(test)]
mod tests {
    use bio::data_structures::annot_map::AnnotMap;
    use crate::reads::as_annotmap::IntoAnnotMap;
    use rust_htslib::bam;
    use std::path::Path;
    use rust_htslib::bam::Read;
    use rust_htslib::bam::Records;
    use bio_types::annot::contig::Contig;
    use bio_types::strand::Strand;
    use crate::locus::from_read::*;
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
            .add_to(&mut map, &hd);

        assert_eq!(res, Ok(()))
    }


}
