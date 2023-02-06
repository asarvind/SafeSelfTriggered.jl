# SafeSelfTriggered

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://htmlpreview.github.io/?https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.html)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://htmlpreview.github.io/?https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.html)

This package implements a self triggered controller synthesis algorithm for safe execution from the paper titled **Safe Self-Triggered Control Based on Precomputed Reachability Sequences**.  This paper will appear in the *26th ACM International Conference on Hybrid Systems: Computation and Control*.  A preliminary version of the paper is available as the file `paper.pdf` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs).

## Documentation
The documentation of the tutorial is available as a Jupyter notebook `tutorial.ipynb` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs) and also this [webpage](https://htmlpreview.github.io/?https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.html).  Clone the repository into the host machine.  Then execute the cells in the tutorial notebook `docs/tutorial.ipynb` or code blocks in this [webpage](https://htmlpreview.github.io/?https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.html), to understand the methods in the package.

## Installation 

### Mosek License
You will need a Mosek license for executing the programs in this repo.  Download a Mosek license and include it in a folder `mosek` in the home directory, like `$Home/mosek/mosek.lic`.  Detailed instructions are provided in this [webpage](https://docs.mosek.com/latest/install/installation.html).  

### Precompiling
Clone the repository anywhere on the host machine. Then from any working folder, do the following to access the methods of the package from that working folder.
```julia   
using Pkg;
Pkg.activate("full_path_to_the_repository_folder"); # Eg. Pkg.activate("/Users/arvind/main/programs/devpackages/SafeSelfTriggered.jl")
Pkg.instantiate()
Pkg.resolve()
```

*Note: This package is not yet registered in the official Julia registry.  So, you can not install it through the Julia registry.*