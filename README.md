# MAPiT
## MAPiT: measure-preserving *MAP* of pseudotime *i*nto true *T*ime
This repository contains a Matlab implementation of MAPiT for
the transformation of pseudotime data to meaningful scales. Application cases in
the associated manusucript "*Reconstructing temporal and spatial dynamics in single-cell experiments*" are
demonstrated in example code.

![MAPiT](/figs/MAPiT.png)    

## Feature overview
MAPiT features include

* Transformation of pseudotime trajectories to real-time scales
* Reconstruction of spatial arrangement of tumor spheroids 
* Easy integration with pseudotime algorithms


## Installation
MAPiT itself is not a software package that has to be installed, but consists of a set of Matlab scripts that have to be adapted and called within a Matlab session for each application case.

MAPiT has no dependencies to other third-party Matlab interfaces/toolboxes.
However, the examples use Wanderlust (part of the Cyt3 toolbox https://github.com/dpeerlab/cyt3) and Diffusion Maps (https://www.helmholtz-muenchen.de/icb/research/groups/machine-learning/projects/dpt/index.html) to derive pseudotime values from single-cell experimental data. Furthermore, kernel density estmation with linked boundary conditions (https://github.com/MColbrook/Kernel-Density-Estimation-with-Linked-BCs.git) is used to obtain distribution in pseudotime for the cell cycle example.

## Examples

This repository contains the scripts necessary to perform the nonlinear
transformation from pseudotime to real-time, or spatial scale with MAPiT, and 
associated workflows for all examples presented in the manuscript. 
The examples are:

- Toy example
[`example_toy.m`](/examples/example_toy.m)
- Cell cycle analysis
[`example_cellcycle.m`](/examples/example_cellcycle.m)
- Spheroid analysis
[`example_spheroid.m`](/examples/example_spheroid.m)
- Paper Fig 2: Cell cycle - MAPiT vs. microscopy
[`MAPiT_CellCycle_Validation.m`](/examples/MAPiT_CellCycle_Validation.m)
- Paper Fig 4: Spheroid - MAPiT vs. microscopy
[`MAPiT_Spheroid_Validation.m`](/examples/MAPiT_Spheroid_Validation.m)

## Usage
Workflow for analysing single-cell data with MAPiT

1. Generate pseudotemporal ordering of cells with your favorite algorithm
2. Define true-scale distribution 
3. Get joint distribution of pseudotime and markers with
[`jointDensityPseudotimeY.m`](jointDensityPseudotimeY.m)
4. Get transformation with [`preMAPiT.m`](preMAPiT.m)
5. Transform pseudotime trajectories to new scale with [`MAPiT.m`](MAPiT.m)

Usage of MAPiT with R or Python based pseudotime analysis methods is straight
forward. Entry point of pseudotime values and single cell data is step 3. 
The function [`jointDensityPseudotimeY.m`](jointDensityPseudotimeY.m) takes a `nx1` vector of
pseudotime values from the pseudotime algorithm and a `nx1` vector of a marker signals from `n` single cells of the single cell
dataset as input. These vectors, must be imported to Matlab, e.g. by importing
 a `.csv` file.

## Citation
**Citeable DOI for the latest MAPiT release:**
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3630379.svg)](https://doi.org/10.5281/zenodo.3630379)

When using MAPiT in your project, please cite

*Reconstructing temporal and spatial dynamics from single-cell pseudotime using prior knowledge of real scale cell densities*
Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm, Frank Allgöwer  
Sci Rep 10, 3619 (2020). https://doi.org/10.1038/s41598-020-60400-z
[![DOI:10.1038/s41598-020-60400-z](https://zenodo.org/badge/DOI/10.5281/zenodo.3630379.svg)](https://doi.org/10.1038/s41598-020-60400-z)
 
