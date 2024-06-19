# Computing interval replacement of persistence modules

We provide a code implementation for computing interval replacements of persistence modules **[1]**. Interval replacements are computed using Möbius inversion of interval ranks, also called compression multiplicities **[2]**. Interval ranks are computed with the formula provided in **[1]** via rank computation of linear maps. 

## Background

In topological data analysis, the usual pipeline is the following:
1. Geometric data
2. Filtered simplicial complexes (obtained with a choice of filtration)
3. Persistence module M (obtained with a choice of homology functor)
4. **Algebraic signature of M (obtained with a choice of algebraic descriptor)**
5. Vector that can be used in machine learning tasks (obtained with a choice of vectorization)

Tasks 1 to 3 should be straightforward to implement thanks to existing code on other repositories. This repository focuses on **task 4**. Given a persistence module, the code computes its interval replacement, used as an algebraic descriptor. This choice of descriptor is motivated by **[1]**: the interval replacement preserves the so-called interval ranks, which generalize the rank invariant property of other algebraic descriptors. Task 5 is still an open problem for the interval replacement.

## Overview

In this repository, you will find:
- a **tutorial** notebook, in which we explain how to compute interval replacements, interval ranks, and other features ;
- a **utils.py** file containing all the necessary code implementations used in the notebook ;
-  and also a **display.py** file containing some useful functions to visualize representations in 1D or 2D settings.

## Key features

In the notebook, you will find how to:
- **instantiate a d-dimensional grid** which is the quiver considered here. This is done within the class `Representation` ; 
- **define a representation (persistence module)** by adding vector spaces and linear maps to the quiver ;
- **define intervals** of the quiver with the class `Interval`. By default, intervals are defined by a list of sources and a list of sinks. One can access all points within the interval using `int_hull`. Conversely, given a list of points forming a connected and convex set, one can instantiate an `Interval` object by using `get_src_snk` ;
- **obtain the list of all intervals** thanks to `list_int` ;
- **compute the interval rank** of a given interval. It is computed with the formula from **[1]** via rank computation of linear maps ;
> [!NOTE]
> Rank computations with `np.linalg` are done in $\mathbb{R}$. You can implement another rank computation algorithm like Gaussian eliminations to implement rank computation in finite fields for example. 
- **compute the interval signed multiplicity** of a given interval via Möbius inversion, by computating the cover of the interval. Signed multiplicities yield the interval replacement of the persistence module.

Additionally, we provide some **visualization** features for the quiver and its intervals.

By default "total (tot)" compression system is used. To use "source-sink (ss)" insted, add a `compression` argument: `L.int_replacement(interval, compression='ss')`.

## Installation

This implementation is built from scratch and does not depend on any external Python libraries, except NumPy and IPython for visualization purpose in the tutorial. You can install them with `pip` by running the following commands:
```
pip install numpy
```
and 
```
pip install ipython
```

## Usage

This code is distributed in the hope that it will be useful. It might be integrated into a topological data analysis pipeline to provide algebraic descriptors from data directly. 

## Future features

- [ ] Implementation of rank computation in finite fields

## References

**[1]**: Asashiba, H., Gauthier, E., & Liu, E. _Interval Replacements of Persistence Modules._ arXiv preprint [arXiv:2403.08308](https://arxiv.org/abs/2403.08308) (2024). 

**[2]**: Asashiba, H., Escolar, E. G., Nakashima, K., & Yoshiwaki, M. _On Approximation of 2D Persistence Modules by Interval-decomposables._ Journal of Computational Algebra, Volumes 6–7, 2023, 100007, ISSN 2772-8277, [https://doi.org/10.1016/j.jaca.2023.100007](https://doi.org/10.1016/j.jaca.2023.100007).

**[3]**: Kim, W., & Mémoli, F. *Generalized persistence diagrams for persistence modules over posets*. J Appl. and Comput. Topology 5, 533–581 (2021). [https://doi.org/10.1007/s41468-021-00075-1](https://doi.org/10.1007/s41468-021-00075-1).

**[4]**: Asashiba, H., Buchet, M., Escolar, E. G., Nakashima, K., & Yoshiwaki, M. *On interval decomposability of 2D persistence modules*, Computational Geometry, Volumes 105–106, 2022, 101879, ISSN 0925-7721, [https://doi.org/10.1016/j.comgeo.2022.101879](https://doi.org/10.1016/j.comgeo.2022.101879).
