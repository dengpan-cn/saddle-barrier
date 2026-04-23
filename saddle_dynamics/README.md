# contributor

This project was developed by members in Lei Zhang's group at Peking University. Please contact Shuonan Wu, Yuchen Xie or Lei Zhang (<zhangl@math.pku.edu.cn>) for further information.

# jampel

This is a MATLAB project for exploring the energy-landscape structure of jammed particle packings with high-index saddle dynamics (HiSD). The current workflow starts from a generated index-4 saddle, searches downward through connected saddles of lower index, and analyzes the resulting network of minima.

## What This Repository Does

The codebase has two main stages:

1. Generate or load saddle configurations in a hierarchy from index 4 down to index 0.
2. Analyze the resulting minima network using geometric distance and barrier-based graph quantities.

In the current version:

- [main/Initialization.m](main/Initialization.m) samples random initial packings, searches for a high-energy minimum, then climbs to index-1 through index-4 saddles.
- [main/main.m](main/main.m) starts from an index-4 saddle and performs repeated downward HiSD searches to recover connected index-3, index-2, index-1, and index-0 states.
- [analysis/jampel_calcsubD.m](analysis/jampel_calcsubD.m) computes pairwise distances between minima and constructs `d_mat` and `sub_d_mat`.
- [analysis/jampel_calcsubE.m](analysis/jampel_calcsubE.m) computes barrier-based connectivity between minima and constructs `sub_E_mat`.

## Repository Structure

```text
jampel/
|-- main/
|   |-- Initialization.m
|   |-- main.m
|   |-- hisd_sirqit.m
|   |-- hisdoptions.m
|   |-- up_search.m
|   |-- down_search_basic.m
|   |-- jampel_calcForce.m
|   |-- jampel_calcHessian.m
|   |-- jampel_calcDistance.m
|   |-- jampel_compareConf.m
|   |-- jampel_cvtRun2ArXiv.m
|   `-- Check_struct.m
|-- analysis/
|   |-- jampel_calcsubD.m
|   `-- jampel_calcsubE.m
`-- README.md
```

## Requirements

The code is written for MATLAB.

Recommended MATLAB functionality:

- `graph`, `minspantree`, `shortestpath`, and `distances`
- `eigs`
- `matchpairs`
- `parfor` if you want parallel acceleration in `analysis/jampel_calcsubD.m`

Depending on your MATLAB installation, some of these functions may require Optimization Toolbox, Parallel Computing Toolbox, or a sufficiently recent MATLAB release.

## Data Flow

The scripts exchange data through `.mat` files stored in `./results`.

Main generated files include:

- `saddle_order_0.mat` to `saddle_order_4.mat`
- `d_mat.mat`
- `sub_d_mat.mat`
- `sub_d_hub.mat`
- `sub_E_mat.mat`

## Suggested Workflow

Run the scripts in the following order from MATLAB:

1. Change into the `main` directory and run `Initialization.m`.
2. Still in `main`, run `main.m`.
3. Change into the `analysis` directory and run `jampel_calcsubD.m`.
4. Run `jampel_calcsubE.m`.

In MATLAB, a typical sequence is:

```matlab
cd main
Initialization
main

cd ../analysis
jampel_calcsubD
jampel_calcsubE
```

## Main Components

### Saddle search

- `hisd_sirqit.m` is the core HiSD solver.
- `hisdoptions.m` creates the solver option structure.
- `up_search.m` perturbs a configuration toward a higher-index saddle.
- `down_search_basic.m` perturbs a configuration along a selected unstable mode and searches downward.

### Configuration utilities

- `jampel_calcForce.m` computes forces, energy, contact number, and pressure-related quantities.
- `jampel_calcHessian.m` computes the Hessian matrix.
- `jampel_compareConf.m` checks whether two archived configurations are equivalent.
- `jampel_calcDistance.m` defines a contact-based distance between two configurations.
- `jampel_cvtRun2ArXiv.m` converts a runtime state into a compact archived structure.
- `Check_struct.m` tests whether a configuration is already present in a set.

### Analysis

- `jampel_calcsubD.m` constructs a graph over minima using pairwise configuration distance.
- `jampel_calcsubE.m` constructs a barrier graph from index-1 saddles and evaluates a reduced `sub_E_mat` on the largest connected component.

## Notes

- To compute minima or saddle-point configurations at different packing fractions or with different interaction potentials, the corresponding parts of the code need to be modified accordingly.

