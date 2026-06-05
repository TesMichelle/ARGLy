# ARGLy Simulation & Validation Scripts

This directory contains Python scripts for validating the MLE engine through simulation and bootstrap analysis.

## Dependencies

Install the required Python packages:
```sh
pip install tskit msprime scipy pandas numpy
```

Additionally, you need SINGER and POLEGON compiled and available. By default, the scripts look for `singer_native` and `polegon_native` inside your system `PATH` or parent research folder. You can override their executable locations using environment variables:
- `SINGER_PATH` (defaults to SINGER C++ binary)
- `POLEGON_PATH` (defaults to POLEGON C++ binary)
- `CONVERT_TSKIT_PATH` (defaults to SINGER's convert_to_tskit script)

## Usage

All scripts should be executed from the repository root directory.

### Run demographic scenarios benchmark:
```sh
python scripts/pipeline.py --scenarios
```

### Run parametric bootstrap (MCMC reconstructed trees):
```sh
python scripts/run_bootstrap.py
```

### Run Oracle bootstrap (True trees):
```sh
python scripts/run_bootstrap_true_trees.py
```
