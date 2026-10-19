# ARGLy
Maximum Likelihood method for detecting ghost population with ARG.

## Build

Preferred simple build: GNU Make, C/C++ compilers, and NLopt with its C++ header
(`nlopt.hpp`). On macOS: `brew install nlopt`. On Debian/Ubuntu:
`apt install libnlopt-dev libnlopt-cxx-dev`. The fit modes use NLopt BOBYQA.

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
proportion, and the base, ghost, and outgroup population sizes. Without
`--grid-*` options, ratio mode evaluates just these fixed values.

To evaluate a grid, add any `--grid-*` axis, for example:

```sh
./build/argly --mode ratio --input example.arg \
  --sample-population 0,0,0,1 --sample-id 0 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --grid-time 900,100,9 --grid-prob 0,1,11 \
  --grid-N-base 10000,20000,5 --grid-N-ghost 12000,24000,5 \
  --output example.ratio-grid.tsv
```

`--grid-time` runs from high to low; other grid axes can run in either
direction. Each population-size axis has either one value or the same number
of values as the other varying size axes. Missing axes use `--params`.
The fast grid pass shares node likelihood calculations across points and
accumulates base and ghost scenarios separately in log space.

The columns are `tree_index`, `left`, `right`, `log_base`, `log_ghost`, and
`log_ratio`. The first two likelihood columns sum on the probability scale;
`log_ratio = log_base - log_ghost` (base divided by ghost). `left` and `right`
identify the genomic interval of each tree. If `--output` is omitted for a
single input file, the TSV is printed to standard output. For batch files,
use `--path DIR --prefix PREFIX --sim-count N [--suffix .arg]`; each input
gets its own `.ratio.tsv` file.
With a grid, there is one TSV row per tree and grid point. The row includes
`time_index`, `N_index`, `prob_index`, the corresponding parameter values,
and the same three log likelihood columns.

## Repeated likelihood queries with a tree ensemble cache

`--mode cached` prepares ALL local trees of one tree sequence and caches
sufficient statistics at their event nodes (plus split boundaries and zero).
Each record contains coalescence counts, pair-time exposures and lineage
counts for the base and ghost populations. Equal-time nodes share a record.
The evaluator combines a structured prefix with the common all-base suffix.
For each population the log likelihood is
`coalescences * log(1/(2*N)) - exposure/(2*N)`, plus the migration weights.
Changing `N` therefore does not rebuild node statistics.

First inspect memory without constructing the scenario cache:

```sh
./build/argly --mode cached --input example.arg \
  --sample-population 0,0,0,1 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --estimate-only
```

`--estimate-only` prints the estimate without allocating scenario records,
even if a saved cache exists. Normal runs report whether statistics were built
or loaded from disk. The estimate reports trees, summed scenario counts, record counts,
statistics bytes and estimated resident cache bytes. It excludes the loaded
tskit tree sequence, allocator overhead and temporary construction buffers;
it is not a prediction of total process peak RSS. Metadata is prepared during
estimation. Statistics storage is approximately
`sum(boundaries_per_tree * 2^free_lineages_per_tree) * sizeof(LikelihoodStatistics)`
(currently 32 bytes per record on the supported build). The cache uses local
tree boundaries; it does not allocate the whole tree-sequence node array per
scenario. `--cache-limit-mb 4096` rejects a cache estimated above 4096 MiB
before its statistics are built or loaded.

To evaluate several parameter sets in one process, put nine numbers per line
in `queries.txt`, comma- or whitespace-separated:

```text
500,1000,2000,0.2,0.5,0.5,10000,12000,9000
501,1000,2000,0.21,0.5,0.5,11000,13000,9000
850,1000,2000,0.15,0.5,0.5,9000,14000,9000
```

```sh
./build/argly --mode cached --input example.arg \
  --sample-population 0,0,0,1 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --queries queries.txt --cache-limit-mb 4096 --output cached.tsv
```

`--params` defines the fixed split times used to build the cache. Every query
may change migration time, migration probability and population sizes, but
must keep both split times unchanged. Times must satisfy
`0 <= migration_time < split2 < split1`. The current model requires binary
trees with one root and present-day sample leaves. As in the existing
likelihood calculation, split proportions and the isolated outgroup size
do not contribute to the score; outgroup-only coalescences below split2 are
omitted. Impossible admixed/outgroup coalescences give `-inf`.

Without `--queries`, the program evaluates `--params` once. Use `--queries -`
to keep the process alive and read successive parameter sets from stdin;
each result is flushed immediately. Blank lines and `#` comments are allowed.
By default, cached mode keeps the statistics in RAM for that process. Add
`--tree-cache FILE` to write all per-tree event boundaries, scenario statistics
and root statistics to disk. A later invocation with the same `FILE` loads
them and can evaluate **new** parameter sets without rebuilding the trees.
The file is binary and tied to its format/ABI. The program rejects it if the
tree-sequence contents, sample populations, population indices or fixed split
parameters have changed. `--estimate-only` still computes the memory estimate
without loading or building the full scenario statistics.

Default TSV output is `query_index, log_llh`, where `log_llh` is the unweighted
sum over ALL local trees. `--per-tree` instead writes
`query_index, tree_index, left, right, log_llh` for each query/tree pair.
`--tree-num` does not subsample cached mode. Preparation, cache construction
and individual query times are reported separately on stderr. When comparing
against grid mode, make sure the same local trees are evaluated.

Each tree remembers its last migration-time interval. New queries first test
that interval and its two neighbours, falling back to binary search for
larger jumps. Queries still sum over all scenarios; cache memory and execution
time grow exponentially with the number of free lineages at split1.

## Finding the maximum

For a separate estimator with piecewise `Nghost(t)` and optional piecewise `Nbase(t)`, see
[the independent epoch ECME module](docs/epoch_fit.md). Build it with
`make build/argly-epoch-fit`; its CLI and caches are separate from the modes below.

The three fit modes maximize the **unweighted sum of log likelihoods across
all trees** in a single input tree sequence. They build the same scenario
cache as `--mode cached`, so use `--estimate-only` or `--cache-limit-mb` when
checking memory. The optimized parameters are migration time, migration
probability, `Nbase`, and `Nghost`. Both split times and `Noutgroup` come from
`--params` and remain fixed. Split proportions and `Noutgroup` currently do
not affect this likelihood, so optimizing them would be meaningless.

Each `--fit-*` axis is `LOWER,UPPER,COUNT`. Time and probability are spaced
linearly; population sizes are spaced logarithmically. Bounds must satisfy
`0 <= time < split2`, `0 <= probability <= 1`, and positive population sizes.
The endpoints are included. Set equal lower and upper bounds, with count 1,
to fix a parameter.

```sh
./build/argly --mode fit-hybrid --input example.arg \
  --sample-population 0,0,0,1 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --fit-time 100,900,9 --fit-prob 0,1,11 \
  --fit-N-base 10000,10000,1 --fit-N-ghost 5000,20000,5 \
  --fit-starts 3 --fit-maxeval 250 --output fit.tsv
```

Choose `fit-grid` for the full Cartesian grid, `fit-opt` for a local NLopt
search starting at `--params`, or `fit-hybrid` for the grid followed by local
searches from up to `--fit-starts` seeds: high-scoring separated grid points
and the supplied `--params` point. Local search
uses BOBYQA in normalized coordinates (N is log-transformed); if only one
parameter varies it uses Nelder–Mead. `--fit-maxeval` applies to each local
start. The TSV records every evaluation and a final `best` row. The same
best value, full nine-parameter vector, number of evaluations, and cache
estimate are printed to stderr. The full grid is capped at one million points.
Add `--tree-num 1000` to fit modes to sample trees at the centers of 1000
equally spaced genomic intervals in the full ARG. If multiple positions fall
in one long tree, that tree is included once, so the number of distinct trees
can be below 1000. If the ARG contains fewer than 1000 trees, all are used.
Without this option, fit modes use every tree. The persistent cache signature
includes the requested sample count and selection method.
Stderr also reports cache preparation, build, load and save progress in percent,
then the estimated cache memory and peak process RSS. Fit runs report the grid
optimum before local search and progress for each optimizer start. TSV output
remains separate from these monitoring messages.

Fit modes automatically save the **tree statistics**, not the scores of grid
points, to `INPUT.tree-cache.SIGNATURE.bin`. The next run with the same input
and fixed parameters loads that file; it still evaluates the likelihood at
every requested grid or optimizer point, including points that were never on
the old grid. The signature depends on the input file contents, sample
populations, population indices, split times/proportions and fixed outgroup
size. Changing time, probability, `Nbase` or `Nghost` ranges reuses the same
file. Use `--tree-cache FILE` to choose another location. A mismatched custom
cache file is rejected. The file is a local build artifact; do not move it
between machines with different binary ABIs.

The C++ API in `likelihood_cache.h` separates estimation, construction and
queries:

```cpp
treellh::TreeEnsemble ensemble(ts, parameters, sample_population, 0, 1);
auto estimate = ensemble.memory_estimate();
ensemble.build();
double total = ensemble.evaluate(new_parameters);
auto per_tree = ensemble.evaluate_trees(another_parameter_set);
```

The tree sequence must remain alive and unchanged while the ensemble exists.
Queries update per-tree cursors, so concurrent queries require separate
ensemble objects or external synchronization. Regression checks are available
with `make test-cache` or `meson test -C builddir`.
