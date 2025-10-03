use criterion::{Criterion, criterion_group, criterion_main};
use std::{hint::black_box, sync::LazyLock};

const SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| {
    std::env::var("SEQ")
        .expect("bench requires a sequence to run")
        .into_bytes()
});

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("dg", |b| {
        b.iter(|| seqfold_rs::fold::dg(black_box(&SEQ), black_box(Some(37.0))))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
