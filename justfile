SEQ := "GAAATAGACGCCAAGTTCAATCCGTACTCCGACGTACGATGGAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGGGTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTT"

perf-core:
    cargo build -r -p seqfold && \
    perf record --call-graph dwarf \
        ./target/release/seqfold {{SEQ}}

time-core:
    cargo build -r -p seqfold && \
    /usr/bin/time -v \
        ./target/release/seqfold {{SEQ}}

