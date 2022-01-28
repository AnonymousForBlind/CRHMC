# Constrained Riemannian Hamiltonian Monte Carlo (CRHMC)

This is a `Matlab` and `C++` implementation of constrained Riemannian Hamiltonian Monte Carlo for sampling from high dimensional disributions on polytopes. This document explains how one can reproduce our experimental results in the submission.

## Folder/File Structure

1. `Instances`: this folder includes constrained-based models from systems biology and LP test sets in `0raw` folder. Running `preparePolytope.m` prepares instances for each algorithm -- it preprocesses constrained-based models by a presolver in the CRHMC package, transforms them into full-dimensional instances, and then round the full-dimensional ones by the maximum volume ellipsoid (MVE) algorithm.
2. `RHMC-test`: `testpaper_rhmc.m` lets CRHMC run on the instances and stores result mat files in `rhmc_test` subfolder.
3. `CHRR-test`: this folder contains CHRR files downloaded from Bounciness/Volume-and-Sampling repository. `testpaper_chrr.m` lets CHRR run on the rounded instances and stores result mat files in `chrr_test` subfolder.
4. `Volesti-test`: this folder contains volEsti-related files. `VolestiSample.R` runs on the rounded instances and stores samples drawn. Then `testpaper_volesti.m` measures quantities of interest based on these samples and stores result mat files in `volesti_test` subfolder.
5. `PolytopeSamplerMatlab-develop`: this folder contains the implementation of CRHMC.
6. `Benchmark`: `testpaper_benchmark.m` runs CRHMC on structured instances such as hypercubes, simplices, and Birkhoff polytopes, and then stores result mat files in subfolders.
7. `testpaper_plot.m`: based on sampling results from the bio and LP dataset, it draws plots for sampling time and mixing rate versus dimension and the number of nonzeros.
8. `testpaper_bench_plot.m`: based on sampling results from the benchmark instances, it draws plots for quantities of interest.

## How to Reproduce Results

We outline the whole steps first and then go over it one-by-one.
1. Data preparation from Bio and LP instances.
2. Run CRHMC, CHRR, and CDHR on these instances.
3. Run CRHMC on the benchmark instances.

### Data Preparation

You should first create two empty subfolders `1chrr` and `2cdhr` under `Instances` folder. Then simply run `preparePolytope.m` to prepare rounded instances for CHRR and CDHR.

### Algorithms on Bio and LP Instances

When running experiment on these instances, execute `MATLAB` by `startMatlab.bat`, which limits all algorithms to one core for fair comparison.

`CRHMC`: we set most of parameters as in `default_options.m` and others as in `testpaper_rhmc.m` (ex. `seed` and `maxTime`). Run `testpaper_rhmc.m`.

`CHRR`: Run `testpaper_chrr.m`.

`CDHR`: After installing `R` and the `volEsti` package, run `VolestiSample.R`. It saves drawn samples in folder `SavedPoints`. Run `testpaper_volesti.m` to obtain mixing rate and sampling time per effective sample, and it saves the measurements for each instance in folder `volesti_test`.

### Algorithms on Benchmark Instances

When running experiment on these instances, execute `MATLAB` by `startMatlab_bench.bat`. Simply run `testpaper_benchmark.m`.
