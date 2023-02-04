# SafeSelfTriggered

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.ipynb)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/devtutorial.ipynb)

This package implements a self triggered controller synthesis algorithm for safe execution from the paper titled **Safe Self-Triggered Control Based on Precomputed Reachability Sequences**.  This paper will appear in the *26th ACM International Conference on Hybrid Systems: Computation and Control*.  A preliminary version of the paper is available as the file `paper.pdf` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs).

The documentation of the tutorial is available as a Jupyter notebook `tutorial.ipynb` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs).  Clone the repository into the host machine.  Then execute the cells in the tutorial notebook `docs/tutorial.ipynb` to understand the methods in the package.

### Installation 
Clone the repository anywhere on the host machine. Then from any working folder, do the following to access the methods of the package from that working folder.
```julia   
using Pkg;
Pkg.activate("full_path_to_the_repository_folder"); 
Pkg.instantiate()
Pkg.resolve()
```

*Note: This package is not yet registered in the official Julia registry.  So, you can not install it through the Julia registry.*