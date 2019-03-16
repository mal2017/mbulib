use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use crate::io::genomic_loci::LocusLike;

// for reference:
// https://stackoverflow.com/questions/30630810/using-generic-iterators-instead-of-specific-list-types
// https://github.com/mal2017/mbulib/commit/690a67d4be92dc3fb3dff8588d955a0e3c7723bb

pub trait IntoAnnotMap : Iterator {
    fn collect_annotmap(&self) -> AnnotMap<String, String> ;
}

impl IntoAnnotMap for Iterator<Item = LocusLike> {
    fn collect_annotmap(&self) -> AnnotMap<String, String>
    {
        let mut map: AnnotMap<String,String> = AnnotMap::new();
        map
    }
}

impl IntoAnnotMap for bam::Records<'_, rust_htslib::bam::Reader> {
    fn collect_annotmap(&self) -> AnnotMap<String, String>
    {
        let mut map: AnnotMap<String,String> = AnnotMap::new();
        map
    }
}


#[cfg(test)]
mod tests {
    use crate::io::genomic_loci::*;
    use crate::io::into_annotmap::IntoAnnotMap;
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
            .records()
            .into_iter()
            .collect_annotmap();


    }
}
