# mtr-phylogeny-mpi

Author: Riki Kawahara<4956309606@edu.k.u-tokyo.ac.jp> and Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>

This repository implements single functionality -- alignments between tandem repeats.

# Requirements

Rust(>1.59)

# Install

After installing rust language, type

```bash
cargo run --release 
```

at the root of this repositoy. Then, `./target/release/align_tandem_repeats` is the binary.

# How to use

We need two files -- reads and units. They should be both (possibly multi-line) fasta files. We assume the reads are HiFi reads from tandem repeats, and the units are the units of the repeats.

`./target/release/align_tandem_repeats --reads <READS> --units <UNITS>` would write the normalized distance matrix to the stdout.
This binary is fully parallelized, and the # of the threads would be automatically estimated.