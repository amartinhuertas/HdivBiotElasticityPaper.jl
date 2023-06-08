# HdivBiotElasticityPaper

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://amartinhuertas.github.io/HdivBiotElasticityPaper.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amartinhuertas.github.io/HdivBiotElasticityPaper.jl/dev/)
[![Build Status](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/amartinhuertas/HdivBiotElasticityPaper.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/amartinhuertas/HdivBiotElasticityPaper.jl)

Repository that holds the Julia software and computational results reported in Section 7.4 and 7.5 of the paper _"Efficient and reliable divergence-conforming methods for an elasticity-poroelasticity interface problem"_ by _S. Badia, M. Hornkjöl, A. Khan, K.-A. Mardal, A. F. Martín, and R. Ruiz-Baier_. While the results in Section 7.1 of the paper were not actually obtained with the software in this repository, the repository also contains Julia scripts to reproduce the results in this section.

In particular:

* The convergence results in Section 7.1 can be generated with the script available [here](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/main/test/ConvergenceTests.jl).
* The T-vessel application problem domain results in Section 7.4 can be generated with the script available [here](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/main/test/TransientHdivBiotElasticityTests.jl)
* The plots in section 7.5 can be generated using the [Pluto.jl](https://github.com/fonsp/Pluto.jl) interactive notebook available [here](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/main/notebooks/HdivBiotElasticityRieszPreconditioningTests.jl). This notebook loads a set of BSON data files available [here](). This data files were in turn automatically generated using the [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) script available [here](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/main/experiments/HdivBiotElasticityRieszPreconditioningTests/run_experiment.jl). These scripts are designed such that the user might seamlessly run new combination of physical and discretization parameter values. See below for instructions. 

## Instructions to generate and visualize Riesz mapping preconditioner evaluation results (Section 7.5)

**IMPORTANT NOTE**: _As a pre-requisite to follow the instructions below, `DrWatson.jl` must be installed in the main julia 
environment, e.g., `v1.8` if you have Julia 1.8. Please install it in the 
main julia environment before following the instructions in the sequel._

1. Select those combinations of parameter-values to be run by editing the dictionary created by the  `generate_param_dicts()` function in the script available [here](https://github.com/amartinhuertas/HdivBiotElasticityPaper.jl/blob/main/experiments/HdivBiotElasticityRieszPreconditioningTests/run_experiment.jl).
2. Run the script:  
   ```bash
   cd experiments/HdivBiotElasticityRieszPreconditioningTests/
   julia run_experiment.jl
   ```
   Note that the `--project=XXX` flag is not required in the call to the `julia` command. 
   The script is smart enough  in order to locate the `Project.toml` of the project in an ancestor directory. Upon completion, the results are generated in the `data` directory 
of the project, in particular in the `data/HdivBiotElasticityRieszPreconditioningTests/commit-ID` folder, where `commit-ID` denotes the latest commit in the repository. If you have changes in your local repository, then the BSON data files will be generated at `data/HdivBiotElasticityRieszPreconditioningTests/commit-ID-dirty`. The script is intelligent enough so that if it is interrupted and re-run, it will only run those combinations which are pending.
3. Visualize the results with the Pluto notebook as follows. Go to the root folder of the project and run the following command:
   ```
   julia --project=.
   ``` 
   Then, in the Julia REPL, run:
   ```julia
   import Pluto 
   Pluto.run(notebook="notebooks/HdivBiotElasticityRieszPreconditioningTests.jl") 
   ```
   this will trigger a web browser tab with the contents of the notebook. Finally, in the notebook, select the data directory with the results you would like to visualize using the _"Select commit ID data directory"_ drop-down list.