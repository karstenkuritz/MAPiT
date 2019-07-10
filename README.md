# MAPiT
## MAPiT: measure-preserving *MAP* of pseudotime *i*nto true *T*ime
This repository contains a Matlab implementation of MAPiT for
the transformation of pseudotime data to meaningful scales. Application cases in
the associated manusucript "*Reconstructing temporal and spatial dynamics in single-cell experiments*" are
demonstrated in example code.

![MAPiT](/figs/MAPiT.jpg)    

## Feature overview
MAPiT features include

* Transformation of pseudotime trajectories to real-time scales
* Reconstruction of spatial arrangement of tumor spheroids 
* Easy integration with pseudotime algorithms


## Installation
MAPiT itself is not a software package that has to be installed, but consists of a set of Matlab scripts that have to be adapted and called within a Matlab session for each application case.

MAPiT has no dependencies to other third-party Matlab interfaces/toolboxes.
However, the examples use Wanderlust (part of the Cyt3 toolbox https://github.com/dpeerlab/cyt3) and Diffusion Maps (https://www.helmholtz-muenchen.de/icb/research/groups/machine-learning/projects/dpt/index.html) to derive pseudotime values from single-cell experimental data. 

## Examples

This repository contains the scripts necessary to perform the nonlinear
transformation from pseudotime to real-time, or spatial scale with MAPiT, and 
associated workflows for all examples presented in the manuscript. 
The examples are:

- Toy example
- Cell cycle analysis
- Spheroid analysis

## Usage
Workflow for analysing single-cell data with MAPiT

1. Generate pseudotemporal ordering of cells with your favorite algorithm
2. Define true-scale distribution 
3. Get joint distribution of pseudotime and markers with
[`jointDensityPseudotimeY.m`](jointDensityPseudotimeYpre.m)
4. Get transformation with [`preMAPiT.m`](preMAPiT.m)
5. Transform pseudotime trajectories to new scale with [`MAPiT.m`](MAPiT.m)

## Citation
*Reconstructing temporal and spatial dynamics in single-cell experiments*  
Karsten Kuritz, Daniela Stöhr, Daniela Maichl, Nadine Pollak, Markus Rehm, Frank Allgöwer  
bioRxiv 697151; doi: https://doi.org/10.1101/697151
 
