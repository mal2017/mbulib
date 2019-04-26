
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use criterion::black_box;

extern crate rqmap;

use rust_htslib::bam;
use std::path::Path;
use rqmap::rqmap::*;


fn reads_into_ngslib() {
    let bampath = Path::new("test/hs.pe.test.bam");
    let bam = bam::Reader::from_path(bampath).unwrap();

    // TODO Work on this test
    let _r = RQMap::from_reader(bam, LibraryType::Unstranded, None);

}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("reads_into_ngslib()", |b| b.iter(|| reads_into_ngslib()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
