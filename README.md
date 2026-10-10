# ARGLy
Maximum Likelihood method for detecting ghost population with ARG.

## Build

Preferred simple build: GNU Make and C/C++ compilers.

```sh
make
```

This builds `build/argly` and `build/main`. The Makefile also builds bundled
`tskit` and `kastore` from `src/subprojects/tskit-1.0.0`; if that directory is
missing, it unpacks `src/subprojects/packagecache/tskit-1.0.0.tar.xz`, or
downloads the same archive as a fallback.

Useful options:

```sh
make CXX=g++ CC=gcc
make CXX=g++-14 CC=gcc-14
make CXXSTD=-std=gnu++20
make install PREFIX="$HOME/.local"
make clean
```

The code requires C++20 because it uses `std::span`; it does not require
`std::mdspan`.

The old Meson build is still available if needed. Requirements: Meson, Ninja,
and a C/C++ compiler.

```sh
meson setup builddir src
meson compile -C builddir
```

The main binaries are created in `builddir/`:

```sh
./builddir/argly --help
```

## Likelihoods for one sample

Use `--mode ratio` to calculate the two likelihoods for a specified sample
without running a parameter grid. The command processes every tree in the
input tree sequence and writes one TSV row per tree:

```sh
./build/argly --mode ratio --input example.arg \
  --sample-population 0,0,0,1 --sample-id 0 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --output example.ratio.tsv
```

`--sample-id` is the tskit node ID of the sample whose origin you want to
classify. `--sample-population` lists the observed population ID indexed by
sample node ID. By default, population `0` is admixed and `1` is outgroup;
change these with `--admixed-index` and `--outgroup-index`. `--params` contains
exactly nine fixed values in this order: migration time, second split time,
first split time, migration probability, second split proportion, first split
proportion, and the base, ghost, and outgroup population sizes. No `--grid-*`
options are needed.

The columns are `tree_index`, `left`, `right`, `log_base`, `log_ghost`, and
`log_ratio`. The first two likelihood columns sum on the probability scale;
`log_ratio = log_base - log_ghost` (base divided by ghost). `left` and `right`
identify the genomic interval of each tree. If `--output` is omitted for a
single input file, the TSV is printed to standard output. For batch files,
use `--path DIR --prefix PREFIX --sim-count N [--suffix .arg]`; each input
gets its own `.ratio.tsv` file.
