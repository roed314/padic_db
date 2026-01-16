This repository contains code, scripts, data and examples in support of

* The paper "Families of p-adic fields" by Jordi Guàrdia-Rubies, John W. Jones, Kevin Keating, Sebastian Pauli, David P. Roberts, and David Roe (published in the conference proceedings of Lucant 2025)
* The paper "Distinguished defining polynomials for extensions of p-adic fields" by Jordi Guàrdia-Rubies, John W. Jones, Kevin Keating, Sebastian Pauli, David P. Roberts, and David Roe (submitted to ANTS XVII)
* Revisions to the [LMFDB database of p-adic fields](https://olive.lmfdb.xyz/padicField/)

# Overview

* [family_generation.py](family_generation.py) is a Python file containing functions used for adding family data to the LMFDB.
* [AllExtensions.m](AllExtensions.m) is a Magma file containing intrinsics for constructing all extensions with given invariants (degree, ramification polygon, residual polynomial).
* [GenResPack.m](GenResPack.m) is a Magma file containing intrinsics for computing indices of inseparability and ramification polygons for p-adic extensions, and residual polynomials for the generic polynomial associated to a family of p-adic fields.
* [JumpSetPack.m](JumpSetPack.m) is a Magma file containing intrinsics for computing jump sets of p-adic extensions.
* [getslopes.m](getslopes.m) is a Magma file containing miscellaneous intrinsics used in adding fields to the LMFDB.
* [padic_galois_counts.m](padic_galois_counts.m) is a Magma file containing intrinsics for counting p-adic fields with a given Galois group (assuming p=2 and the group has 2-power order)
* [padic_filtrations.m](padic_filtrations.m) is a Magma file containing intrinsics for exploring different possible ramification filtrations on an abstract group.
* [indices.gp](indices.gp) is a Pari/GP script for computing indices of inseparability.
* The scripts folder contains scripts used to run some of these functions at scale, used for creating upload files to add to the LMFDB.
* The examples folder contains example files illustrating how to use some of the intrinsics defined above.
* The data folder contains some data used by the code above, as well as scripts for downloading necessary data from the LMFDB to make some of the other functions above work.

The `PolRedPadic` folder contains implementations for the algorithms described in the second paper, with the following files:

* [charpoly.m](PolRedPadic/charpoly.m) is a Magma file for computing characteristic polynomials using Newton relations.
* [montes.m](PolRedPadic/montes.m) is a copy of the +Ideals Magma package written by J. Guardia, J. Montes and E. Nart
* [oystein.m](PolRedPadic/oystein.m) is a Magma file for computing Oystein polynomials
* [rampoly.m](PolRedPadic/rampoly.m) is a Magma file for computing ramification polynomials, ramification polygons and residual polynomials
* [distinguished.m](PolRedPadic/distinguished.m) is a Magma file implementing the overall algorithm for choosing distinguished defining polynomials for extensions of p-adic fields

# Installation

Clone this repository using

```
git clone https://github.com/roed314/padic_db.git
```
After starting Magma within the `padic_db` folder, `AttachSpec("spec")` will make all of the intrinsics defined in the Magma files available.

After starting Sage within the `padic_db` folder, `%attach family_generation.py` will make the functions defined in that file available for use.

# Usage

See the individual files for explanations of the intrinsics, and the `examples` folder for example usage.


