# SafeSelfTriggered

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/tutorial.ipynb)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/asarvind/SafeSelfTriggered.jl/blob/main/docs/devtutorial.ipynb)

This package implements the self triggered controller synthesis algorithm for safe execution from the paper titled **Safe Self-Triggered Control Based on PrecomputedReachability Sequences**.  The paper will appear in the *26th ACM International Conference on Hybrid Systems: Computation and Control* conference.  A preliminary version of the paper is given as the file `paper.pdf` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs).

The documentation of the tutorial is available as a Jupyter notebook `tutorial.ipynb` in the [following link](https://github.com/asarvind/SafeSelfTriggered.jl/tree/main/docs).  Clone the repository into your local machine.  Then follow the steps in the tutorial to understand the methods in the package.

After downloading the repository, we can access the methods in the package from any folder by activating it.  To activate the package, type in the julia REPL `using Pkg; Pkg.activate("full_path_to_the_repository_folder")`.  


*Note: This package is not yet registered as in the official Julia registry.  So, you can not install it through the Julia registry.*