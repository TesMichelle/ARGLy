# Piecewise population sizes: independent ECME module

`argly-epoch-fit` is a separate executable and C++ module. It leaves the existing
`Scenario_Computer`, tree cache, grid, ratio and fit implementations unchanged.
It uses existing public tree classification helpers and NLopt, already required
by the standard fit modes.

## Build and run

```sh
make build/argly-epoch-fit
make test-epoch
```

Meson builds the executable too; `meson test -C builddir` runs the new tests.

```sh
./build/argly-epoch-fit --input full.arg \
  --sample-population 0,0,0,1 \
  --params 500,1000,2000,0.2,0.5,0.5,10000,12000,9000 \
  --epochs 5 --tree-num 1000 \
  --time-bounds 100,900 --prob-bounds 0.001,0.999 \
  --N-bounds 100,1000000 --iterations 30 --maxeval 100 \
  --output epoch-fit.tsv
```

The input population labels and parameters must correspond to the actual file.
For the repository's prepared 1000-tree benchmark, the tested command is:

```zsh
./build/argly-epoch-fit \
  --input benchmarks/results/cache_1000/selected_1000.trees \
  --sample-population "$(printf '2,%.0s' {1..199})2" \
  --params 1724.1379310344828,2093.103448275862,19275.862068965518,0.05,0.5,0.5,10000,10000,10000 \
  --admixed-index 2 --outgroup-index 0 \
  --epochs 5 --time-bounds 800,2000 --prob-bounds 0.001,0.2 \
  --N-bounds 1000,100000 --iterations 30 --maxeval 60 \
  --output /private/tmp/epoch-fit.tsv
```

`--params` uses the existing nine-value order. By default `Nbase` stays fixed.
Pass `--fit-base` to estimate both `Nbase(t)` and `Nghost(t)` on the same K equal
epochs covering `[0,split1]`. For example, add `--fit-base` to either command above.
`--base-sizes n0,n1,...` supplies K starting base sizes and also enables fitting
them; without it, every base epoch starts at `Nbase` from `--params`.
`--ghost-sizes n0,n1,...` supplies K initial values, youngest epoch first; otherwise
all epochs start at the supplied constant `Nghost`. Each value must be within
`--N-bounds`. Split times are fixed. Above split1 the common ancestral contribution
uses fixed Nbase from `--params`. Below migration time all lineages use the
estimated base profile when `--fit-base` is enabled.
Pass `--fixed-ghost` to hold every ghost epoch at its supplied value while
fitting the base profile and migration parameters. With no `--ghost-sizes`,
all ghost epochs stay at `Nghost` from `--params`.

`--epoch-boundaries t0,...,tK` replaces the equal epochs with explicit shared
boundaries for both populations. Supply exactly K+1 finite, strictly increasing
times; the first must be 0 and the last must equal `split1` from `--params`.
An internal boundary may coincide with `split2`, even when `split2` is only a
placeholder because the input has no outgroup samples.

If the input tree sequence has no samples labelled with `--outgroup-index`,
`split2` in `--params` is a placeholder: it is excluded from the epoch cache
and does not bound migration time. The allowed range is then `0 <= T_adm < split1`,
including when `--time-bounds` is omitted. With outgroup samples, the usual
`0 <= T_adm < split2` condition remains.

`--tree-num K` selects trees containing the midpoints of K equal genomic windows,
with repeated trees counted once. If K is at least the number of trees, all trees
are used. Omitting this option also selects all trees. Scores are unweighted sums
of local tree log likelihoods, as in the existing fit modes.

`--evaluate-only` evaluates a supplied curve without fitting. `--time-bounds` and
`--prob-bounds` can fix a variable by specifying equal endpoints. Starts must lie
inside the bounds. Probability bounds default to `[1e-6,1-1e-6]` to avoid starting
with exactly zero posterior mass for a class of colourings. Explicit 0/1 endpoints
are supported, but EM can get stuck with zero posterior weights there.

## Two cache levels

1. Each selected tree has a common descending array of boundaries: internal node
   times, split times and equal epoch boundaries. Each colouring stores the
   coalescence counts at the upper boundary and the two lineage counts on the
   interval below it. The current `Step` record occupies 16 bytes. Durations and
   epoch indices are shared across colourings.
2. Each colouring stores one accumulated log likelihood at each boundary, before
   that boundary's coalescence. This costs another 8 bytes per record and is rebuilt
   when either N curve changes. Moving migration time uses the neighbouring
   cached values and the partial interval; it does not rebuild the cache.

Intervals at a vertex/epoch boundary use the younger epoch's N for a coalescence
exactly on that boundary: event epochs have convention `(left_time,right_time]`.
A coalescence exactly at migration time belongs to the common all-base tail,
matching the existing cache. Sample times must be zero and trees must be binary
with a single root. Outgroup lineage handling and omission of its separate
likelihood match the current model.

The cache is in memory for the duration of this process. No persistent cache
format or changes to the existing on-disk cache were introduced. Building scans
the ARG and reports progress; `--cache-limit-mb` (default 1024 MiB, 0 disables it)
checks cache array storage before allocation. It does not cap the input tree
sequence, temporary construction arrays, allocator overhead or total RSS.

## Estimation

At each iteration, the E-step computes a normalized weight for each colouring
from the likelihood of the **whole tree**, with the current curve and migration
parameters. Common root and tail contributions cancel in these weights.

For each population, coalescences and pair-time exposure above migration time
are collected in each epoch and averaged over these weights, then summed over
trees. Base statistics also include the common tree segment below migration
time, where all lineages are base. At fixed migration parameters, the M-step sets

```
Npop[e] = expected_exposure[pop,e] / (2 * expected_coalescences[pop,e])
```

and clamps it to the supplied bounds. A positive exposure with zero expected
coalescences selects the upper bound. If both statistics are zero, the old value
is retained: that epoch contains no information about that population size at
this iteration. Without `--fit-base`, only ghost sizes are updated.

After rebuilding the second cache level, NLopt searches time and probability to
increase the **observed** full-tree likelihood with both N profiles fixed. Because this
step optimizes observed likelihood rather than the EM auxiliary function, the
algorithm is ECME. BOBYQA uses normalized coordinates; when just one migration
parameter is free, Nelder–Mead is used. The best evaluated point is retained,
including the starting point. A material decrease in likelihood raises an error;
tiny roundoff decreases roll back the iteration. `--maxeval` is per migration
search, and `--iterations` limits outer iterations. Convergence uses absolute
log-likelihood improvement (`--tolerance`, default 1e-5).

This is a local method. Repeat with different starting curves and migration
parameters to compare optima. There is no smoothing penalty; many short epochs
can yield weakly identified or boundary estimates. Local trees from one ARG are
correlated: the summed score follows existing behaviour and does not provide
independent-tree uncertainty estimates.

## Outputs and monitoring

`--output FILE` writes an iteration TSV with full log likelihood, time,
probability, improvement and all estimated N values. Without it the TSV goes to stdout.
Stderr reports cache progress, array memory, peak process RSS, each migration
search, outer iterations, final convergence status and elapsed time.

With `--output FILE`, the final epoch table defaults to `FILE.epochs.tsv`.
`--profile-output PATH` overrides this or enables a profile without iteration
output. The profile includes epoch boundaries, Nghost and, when fitted, Nbase,
their expected coalescences, exposure and statuses:

- `below_migration`: this entire epoch is below the final migration time;
- `no_information`: both ghost sufficient statistics are zero;
- `upper_bound` / `lower_bound`: Nghost equals the corresponding fit bound;
- `active`: the epoch has information and its size is inside the bounds.

The first two statuses must not be interpreted as estimated historical sizes.
Expected counts may be fractional because they average over latent colourings.

## Verification

`epoch_test` compares constant curves to the existing cache at probabilities 0,
1 and interior values, several N values and migration times including exact node
and epoch boundaries. Variable curves are checked against an independent
two-tip closed form and numerical score derivatives. Tests also check the
N-update and full ECME monotonicity and the memory limit.
