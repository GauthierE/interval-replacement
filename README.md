# Computing interval replacement of persistence modules

We provide a code implementation for computing interval replacements of persistence modules **[1]**. Interval replacements are computed using Möbius inversion of interval ranks, also called compression multiplicities **[2]**. Interval ranks are computed with the formula provided in **[1]** via rank computation of linear maps. 

## Table of contents

## Background

In topological data analysis, the pipeline is the following:
Geometric data --[Filtration]--> Filtered simplicial complexes --[Homology]--> Persistence module M --[Algebraic descriptor]--> Algebraic signature of M --[Vectorization]--> Vector that can be used in machine learning tasks

## Overview

In this repository, you will find:
- A tutorial notebook, in which we explain how to compute interval replacements, interval ranks, and other features. 
- A **.py** file containing all the necessary code implementations used in the notebook.

## Key features

In the notebook, you will find how to:
- **Instantiate a d-dimensional grid** which is the quiver considered here. This is done within the class `Representation`. 
- **Define a representation (persistence module) M** by adding vector spaces and linear maps to the quiver.
- **Define intervals** of the quiver with the class `Interval`. By default, intervals are defined by a list of sources and a list of sinks. One can access all points within the interval using `int_hull`. Conversely, given a list of points forming a connected and convex set, one can define an `Interval` instance by using `get_src_snk`.
- **Obtain the list of all intervals** thanks to `list_int`.
- **Compute the interval rank** of a given interval. It is computed with the formula from **[1]** via rank computation of linear maps. 
- **Compute the interval signed multiplicity** of a given interval via Möbius inversion, by computating the cover of the interval. Signed multiplicities yield the interval replacement of the persistence module.

Additionally, we provide some **visualization** features for the quiver and its intervals.

# Installation

## Usage

## Future features

- [ ] The interval rank computation relies on rank computation of linear maps. This computation is by default implemented in \mathbbR => implement rank computation in finite fields
- [ ] The interval rank formula holds for compression systems with essential covering property => implement interval rank computation for other compression systems, like the _source-sink_ assignment

## References

**[1]**: Asashiba, H., Gauthier, E., & Liu, E. _Interval Replacements of Persistence Modules._ arXiv preprint [arXiv:2403.08308](https://arxiv.org/abs/2403.08308) (2024). 

**[2]**: Asashiba, H., Escolar, E. G., Nakashima, K., & Yoshiwaki, M. _On Approximation of 2D Persistence Modules by Interval-decomposables._ Journal of Computational Algebra, Volumes 6–7, 2023, 100007, ISSN 2772-8277, [https://doi.org/10.1016/j.jaca.2023.100007](https://doi.org/10.1016/j.jaca.2023.100007) (2023).
