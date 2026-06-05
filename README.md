# ARGLy
Maximum Likelihood method for detecting ghost population with ARG.

## Build

Requirements: Meson, Ninja, and a C/C++ compiler.

```sh
meson setup builddir src
meson compile -C builddir
```

The main binaries are created in `builddir/`:

```sh
./builddir/argly --help
```

## MLE Inference & Validation Pipelines

We provide scripts for demographic simulation, ARG reconstruction validation, and bootstrap analysis in the `scripts/` directory:

- `pipeline.py`: Simulates data under different migration scenarios using `msprime`, runs `SINGER`/`POLEGON`, and optimizes parameters using the C++ MLE engine.
- `run_bootstrap.py`: Performs parametric bootstrap ($B=100$ replicates) to estimate standard errors and $95\%$ confidence intervals.
- `run_bootstrap_true_trees.py`: Runs an Oracle bootstrap using true trees (recombination = 0) to validate the theoretical convergence of the MLE model.

See `scripts/README.md` for dependencies and usage details.
