
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use criterion::black_box;

extern crate rqmap;

use rust_htslib::bam;
use std::path::Path;
use rqmap::rqmap::*;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;

fn rqmap_from_reader() {
    let bampath = Path::new("test/hs.pe.test.bam");
    let bam = bam::Reader::from_path(bampath).unwrap();

    // TODO Work on this test
    let _r = RQMap::from_reader(bam, LibraryType::Unstranded, None, None);

}

fn rqmap_from_indexed() {
    let bampath = Path::new("test/hs.pe.test.bam");
    let bam = bam::IndexedReader::from_path(bampath).unwrap();

    let c0: Contig<String,ReqStrand> = Contig::new("chr1".to_string(),
                                                 564475,
                                                 60,
                                                 ReqStrand::Forward);

    // TODO Work on this test
    let _r = RQMap::from_indexed(bam, vec!(c0), LibraryType::Unstranded, None, None);

}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("rqmap_from_reader()", |b| b.iter(|| rqmap_from_reader()));
    c.bench_function("rqmap_from_indexed()", |b| b.iter(|| rqmap_from_indexed()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
