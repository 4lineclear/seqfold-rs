use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const SEQ: &[u8] = b"GAAATAGACGCCAAGTTCAATCCGTACTCCGACGTACGATGGAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGGGTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTT".as_slice();

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("dg", |b| {
        b.iter(|| seqfold_rs::fold::dg(black_box(SEQ), black_box(Some(37.0))))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
