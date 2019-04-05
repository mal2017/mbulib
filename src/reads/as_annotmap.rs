use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use rust_htslib::bam::Records;
use rust_htslib::bam::HeaderView;
use crate::locus::from_read::*;
use crate::utility::scaffold_dict::ScaffoldDict;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;
use rust_htslib::bam::Read;

// for reference:
// https://stackoverflow.com/questions/30630810/using-generic-iterators-instead-of-specific-list-types

/// Trait for converting bam/sam Records into annotation maps.
pub trait AsAnnotMap {
    fn collect_annotmap(self) -> AnnotMap<String, Contig<String, ReqStrand>> ;
    fn add_to(self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>) -> Result<(), &'static str>;
}

impl<T: Read> AsAnnotMap for T {
    fn collect_annotmap(mut self) -> AnnotMap<String, Contig<String,ReqStrand>> {
        let mut map: AnnotMap<String,Contig<String,ReqStrand>> = AnnotMap::new();
        let hd: &HeaderView = self.header();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);

        let z = self.records()
                    .map(|a| a.unwrap())
                    .map(|a| Contig::from_read(&a, false, &sd));

        for i in z.into_iter() {
            map.insert_loc(i);
        }

        map
    }
    fn add_to(mut self, am: &mut AnnotMap<String, Contig<String,ReqStrand>>) -> Result<(), &'static str> {
        let hd: &HeaderView = self.header();
        let sd: ScaffoldDict = ScaffoldDict::from_header_view(hd);
        let z = self.records()
                    .map(|a| a.unwrap())
                    .map(|a| Contig::from_read(&a, false, &sd));

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
    use crate::reads::as_annotmap::AsAnnotMap;
    use rust_htslib::bam;
    use std::path::Path;
    use rust_htslib::bam::Read;
    use rust_htslib::bam::Records;
    use bio_types::annot::contig::Contig;
    use bio_types::strand::Strand;

    #[test]
    fn bam_reads_2_annotmap() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::Reader::from_path(bampath).unwrap();

        let _res = bam
            .collect_annotmap();

    }

    #[test]
    fn annotmap_functioning() {
        let bampath = Path::new("test/hs.pe.test.bam");
        let bam = bam::Reader::from_path(bampath).unwrap();

        let am = bam
            .collect_annotmap();


        let query = Contig::new("chr1".to_owned(),
                                564477,
                                (565723 - 564477),
                                Strand::Unknown);
        let ov = am.find(&query);

        let n = ov.into_iter().count();

        assert_eq!(n, 9);
    }
}
