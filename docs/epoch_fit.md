# Piecewise population sizes: independent ECME module

`argly-epoch-fit` is a separate executable and C++ module. Its likelihood
cache and ECME loop are separate from the original tree cache and fit modes.
It reuses `Scenario_Computer`'s tree-classification helpers and NLopt. The
`argly --mode ratio` implementation also accepts the fitted N profiles through
its own CLI arguments.

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
  --base-epochs 3 --ghost-epochs 5 --tree-num 1000 \
  --time-bounds 100,900 --prob-bounds 0.001,0.999 \
  --N-bounds 100,1000000 --iterations 30 --restarts 20 --maxeval 100 \
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

`--params` uses the existing nine-value order. Without grid options, the model
estimates one constant Nbase and one constant Nghost on `[0,split1]`. Use
`--base-epochs K` and `--ghost-epochs L` independently to estimate K and L
piecewise-constant sizes. An omitted population count remains one. Both
populations are fitted by default; `--fixed-base` or `--fixed-ghost` keeps the
corresponding supplied profile fixed. The old `--fit-base` flag remains accepted
for compatibility but is no longer needed. `--base-sizes` and `--ghost-sizes`
supply starting values, youngest epoch first, with one value per population
epoch; omitted values start at the corresponding N in `--params`. All starting
sizes must lie within `--N-bounds`. Above split1, the ancestral contribution
still uses fixed Nbase from `--params`; below migration time all lineages use
the base profile. Split times are fixed.

With separate grids, adaptive boundaries are the default whenever a population
has more than one epoch and no explicit boundaries. The default
`--epoch-balance count` counts each coalescence at time `0 < t <= split1`
once in **every local tree in the input**, including outgroup-only events. This
scan runs before the ECME cache is built, even when `--tree-num` later limits
the trees used for likelihood fitting. These
same raw event counts supply quantiles for both grids; different requested
epoch counts can therefore give different cuts. The raw counts do not require
knowing which latent migration colouring occurred. This choice aims for an equal
number of observed coalescences per epoch, independently of the initial `T`,
`p` and `N`. It does not guarantee equal *population-specific expected* counts.

`--epoch-balance expected` retains the alternative posterior-weighted method.
It chooses each population's cuts from its expected coalescence counts under
the initial `T,p` and constant Nbase/Nghost in `--params`. Supplied starting
`--base-sizes` and `--ghost-sizes` apply after the boundaries are chosen.
For either adaptive method, cuts fall strictly between distinct event times.
Tied events cannot be split, so counts can only be approximately equal. More
epochs than distinct informative event times is an error. Boundaries remain
fixed through ECME; the expected counts can change as parameters change.

`--base-boundaries 0,t1,...,split1` and
`--ghost-boundaries 0,u1,...,split1` set the two grids manually. Their lengths
determine the corresponding epoch counts unless `--base-epochs` or
`--ghost-epochs` is also given, in which case the count must match. Values must
be finite and strictly increasing. A manual grid can be combined with an
adaptive grid for the other population. An internal boundary may coincide with
a coalescence or `split2`; a coalescence on that boundary belongs to the
younger interval `(left,right]`, not the older one.

`--epochs K` preserves the shared-grid mode: both populations get K equal-duration
epochs from 0 to split1, as before. `--epoch-boundaries t0,...,tK` replaces those
with K+1 explicit shared boundaries; it requires `--epochs` and cannot be
combined with separate epoch counts or boundaries. `--epoch-balance` applies
only to separate grids. Shared explicit boundaries obey the same endpoint and
ordering rules.

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
   times, split times and the union of both population grids. Each interval
   records the base and ghost epoch indices independently. Each colouring stores the
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

The cache is in memory for the duration of this process. In separate-grid mode
with `--epoch-balance expected`, a one-epoch pilot cache is built to compute
expected coalescence-time masses, then released before building the final cache.
The default raw-count scan does not build a pilot cache. Fully manual grids
skip both boundary-selection passes; the initial window diagnostic still scans
the input trees once.
The number of interval records grows with the **sum** of distinct boundaries
of the two grids, not their Cartesian product. No persistent cache
format or changes to the existing on-disk cache were introduced. Building scans
the ARG and reports progress; `--cache-limit-mb` (default 1024 MiB, 0 disables it)
checks cache array storage before allocation. It does not cap the input tree
sequence, temporary construction arrays, allocator overhead or total RSS.

## Estimation

At each iteration, the E-step computes a normalized weight for each colouring
from the likelihood of the **whole tree**, with the current curve and migration
parameters. Common root and tail contributions cancel in these weights.

For each population, coalescences and pair-time exposure above migration time
are collected in the corresponding population's epochs and averaged over these weights, then summed over
trees. Base statistics also include the common tree segment below migration
time, where all lineages are base. At fixed migration parameters, the M-step sets

```
Npop[e] = expected_exposure[pop,e] / (2 * expected_coalescences[pop,e])
```

and clamps it to the supplied bounds. A positive exposure with zero expected
coalescences selects the upper bound. If both statistics are zero, the old value
is retained: that epoch contains no information about that population size at
this iteration. `--fixed-base` or `--fixed-ghost` skips the corresponding N update.

After rebuilding the second cache level, NLopt searches time and probability to
increase the **observed** full-tree likelihood with both N profiles fixed. Because this
step optimizes observed likelihood rather than the EM auxiliary function, the
algorithm is ECME. BOBYQA uses normalized coordinates; when just one migration
parameter is free, Nelder–Mead is used. Each ECME iteration runs 10 local searches
by default. `--restarts K` changes the number of searches (for example,
`--restarts 20`). Each search starts at an independent pseudorandom point sampled
uniformly within the supplied time and probability bounds; fixed coordinates
remain fixed. The sequence is reproducible with the default `--restart-seed 5489`;
change that seed to draw different starts. The current ECME point is also retained
as a candidate, without consuming one of the random searches. All searches share
the same tree and prefix caches and the same fixed N curves; the best evaluated
point is selected. A material decrease in likelihood raises an error; tiny
roundoff decreases roll back the iteration. `--maxeval` is **per restart**, so the
default search budget is up to 10 × 100 evaluations per ECME iteration.
`--iterations` limits outer iterations. Convergence uses absolute log-likelihood
improvement (`--tolerance`, default 1e-5).

This is still a local method. The restarts improve the conditional search over
time and probability at fixed N, but do not replace full ECME runs with different
starting N curves and migration parameters. There is no smoothing penalty; many short epochs
can yield weakly identified or boundary estimates. Local trees from one ARG are
correlated: the summed score follows existing behaviour and does not provide
independent-tree uncertainty estimates.

## Outputs and monitoring

`--output FILE` writes an iteration TSV with full log likelihood, time,
probability, improvement and all estimated N values. Without it the TSV goes to stdout.
`--kernel-output FILE` additionally writes a posterior-weighted Gaussian
estimate of `Nbase(t)` and `Nghost(t)` **once, after fitting**. It does not feed
the smoothed curve back into ECME. `--kernel-bandwidth H` sets the Gaussian
bandwidth in generations (default `split1/50`); `--kernel-points K` sets the
number of evenly spaced output times from 0 through `split1` (default 201).
The default bandwidth is only a plotting convenience; the shape of the curve
depends on this choice, so report the chosen value when using the estimate.
The final migration time is also inserted as an output point if absent.
The long-format TSV columns are `population,time,N_kernel,kernel_C,kernel_E,bandwidth,status`.
`kernel_C` is the posterior-weighted sum of coalescences multiplied by
`exp(-(event_time-time)^2/(2H^2))`; `kernel_E` applies that same kernel to
pair-time exposure, integrating it over each interval with constant lineage
counts. Their ratio gives `N_kernel = kernel_E/(2*kernel_C)`. The kernel is
unnormalized because its common factor cancels. `status` is `ok`, `low_events`
when `kernel_C<1`, `no_events`, `no_exposure`, or `below_migration`. The last
three statuses have `N_kernel=NA`; ghost points at or below the fitted
migration time are always `below_migration`. Near the time-domain edges the
kernel uses only available exposure and events, so the effective neighbourhood
is asymmetric. The output is a conditional descriptive estimate under the
fitted stepwise model, not a new ECME optimum. `--evaluate-only` also supports
this output using the supplied parameters.
Stderr reports cache progress, array memory, peak process RSS, each migration
search, outer iterations, final convergence status and elapsed time. Before the
first ECME iteration it also prints one diagnostic row per population and epoch:
`raw_C`, `raw_mean_K`, `raw_lineage_time`, `raw_pair_time`, `expected_C`, initial
size `N0`, and the idealized variance and standard deviation of its size
estimate. The `raw_*` fields are calculated once from the uncoloured local
trees across **all input trees**, including outgroup events, before ECME starts.
They are repeated for each population's time grid but are not population
specific. For a window of width `w`, `raw_lineage_time` is the sum of
`∫K(t)dt` over the trees, `raw_mean_K = raw_lineage_time/(tree_count × w)`,
and `raw_pair_time` is the sum of `∫K(t)(K(t)−1)/2 dt`. Thus the latter is
the raw pair-time exposure, not an expected population-specific exposure.
These fields characterize the supplied windows and are never recomputed on
ECME iterations. `raw_C` counts all coalescences in each time window.
`expected_C` is
the E-step expectation for that population across the selected fitting trees
at the initial `T,p,N` values. The report gives the mean raw and expected count
per epoch, plus the minimum and maximum raw counts. Events on a window's upper
boundary are included in that window.

The displayed proxy uses `ideal_var_N = N0^2 / expected_C` and
`ideal_sd_N = N0 / sqrt(expected_C)`; it prints infinity when `expected_C=0`.
For independent, fully assigned coalescent waiting times with a fixed number
of events, this follows from the variance of the sum of exponential waiting
times in `Nhat=E/(2C)`. Here the colouring is latent, epochs may censor waiting
times, and local trees of an ARG are dependent. Thus the reported variance is
an optimistic diagnostic, **not** an uncertainty estimate or confidence
interval for the ECME result.

With `--output FILE`, the final epoch table defaults to `FILE.epochs.tsv`.
`--profile-output PATH` overrides this or enables a profile without iteration
output. Shared-grid mode retains the original wide rows with both populations.
Separate-grid mode uses one row per population and epoch with columns
`population,epoch,left_time,right_time,N,expected_coalescences,exposure,status`.
To classify a sample with the fitted sizes, pass each population's
`left_time`/`right_time` boundaries and `N` values to `argly --mode ratio` as
`--base-boundaries`/`--base-sizes` and `--ghost-boundaries`/`--ghost-sizes`.
Keep the other model parameters consistent with the fit.
Statuses include:

- `below_migration`: this entire epoch is below the final migration time;
- `no_information`: both sufficient statistics for that population and epoch are zero;
- `fixed`: the population profile was held fixed;
- `upper_bound` / `lower_bound`: fitted N equals the corresponding fit bound;
- `active`: the epoch has information and its size is inside the bounds.

The first two statuses must not be interpreted as estimated historical sizes.
Expected counts may be fractional because they average over latent colourings.

## Verification

`epoch_test` compares constant curves to the existing cache at probabilities 0,
1 and interior values, several N values and migration times including exact node
and epoch boundaries. Variable curves are checked against an independent
two-tip closed form and numerical score derivatives. Distinct base/ghost grids
are checked against the legacy constant-size likelihood, their expected event
masses against E-step counts, and their score derivatives numerically. Tests also
check balanced boundaries, N-updates, full ECME monotonicity and the memory limit.

`Ensemble::derivatives()` and `derivatives_by_tree()` compute analytic
derivatives of the observed log likelihood with respect to migration time and
probability from the existing scenario and prefix cache. They require
`0<p<1` and a migration time strictly inside a cache interval; the time
derivative is not returned at an event or epoch boundary. The optimizer does
not use these APIs. Synthetic tests compare them with central differences,
Local plotting scripts under `docs/figures` show individual test trees and all
1000 trees of the benchmark input; that generated figure directory is ignored
by Git and is not part of this commit.
