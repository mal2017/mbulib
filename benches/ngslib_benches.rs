
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use criterion::black_box;

extern crate ngslib;

use bio::data_structures::annot_map::AnnotMap;
use rust_htslib::bam;
use std::path::Path;
use rust_htslib::bam::Read;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use ngslib::ngslibrary::*;
use ngslib::locus::shift::*;
use ngslib::utility::scaffold_dict::*;
use rust_htslib::bam::HeaderView;
use bio_types::strand::ReqStrand;


fn reads_into_ngslib() {
    let bampath = Path::new("test/hs.pe.test.bam");
    let bam = bam::Reader::from_path(bampath).unwrap();

    // TODO Work on this test
    let _r = NGSLibrary::from_reader(bam, LibraryType::Unstranded, None);

}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("reads_into_ngslib()", |b| b.iter(|| reads_into_ngslib()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
