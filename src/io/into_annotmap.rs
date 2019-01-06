use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;

// for reference:
// https://stackoverflow.com/questions/30630810/using-generic-iterators-instead-of-specific-list-types
// https://github.com/mal2017/mbulib/commit/690a67d4be92dc3fb3dff8588d955a0e3c7723bb

pub trait IntoAnnotMap {
    fn collect_annotmap(&self) -> (); //AnnotMap<String, String> ;
}

impl IntoAnnotMap for Iterator<Item = bam::Record> {
    fn collect_annotmap(&self) -> () {
        ()
    }
}
