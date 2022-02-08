# RSBilinear
Data from Numerical Trials for Forthcoming Manuscript: "Improved convex and concave relaxations of composite bilinear forms" (hyperlink and citation information will be added once the paper becomes available online)

# Authors
[Matthew E. Wilhelm](https://scholar.google.com/citations?user=sr4baQ0AAAAJ&hl=en&authuser=1), [Matthew D. Stuber](https://cbe.engr.uconn.edu/person/matthew-stuber/)

# Overview of Paper
The generalized McCormick relaxations allows for convex/concave relaxations of nonlinear functions to be computed in the original problem space original problem (without a need to introduce auxillary variables introducing auxiliary variables). We have adapted a [recent contribution](https://link.springer.com/article/10.1007/s10107-020-01541-x) illustrating how additional nontrivial inequality constraints may be used in factorable programming to tighten relaxations of the ubiquitous bilinear term to develop a novel approach for computiong composite relaxations of the bilinear term in the reduced-space. We then detailed three different approaches to generating the necessary apriori relaxations needed to implement this approach: using of a McCormick relaxations coupled to affine arithmetic, the propagation of affine relaxations (derived from subgradients), and direct enumeration of each factor. We include two case studies highlighting the improvements resulting from this approach: a supply-chain based response surface model (RSM) and a dynamic parameter estimation problem. A randomly generated benchmarking method is then used to compare these new approaches to state-of-the-art nonconvex optimizers.

# How to Reproduce Benchmark Results
A CSV file containing the results may be reproduced using in the following fashion. 
1. Clone this git repository to you local machine.
2. Install Julia 1.6 or greater, see [https://julialang.org/downloads/](https://julialang.org/downloads/).
3. Install SCIP and SCIP.jl per the instructions provided at [https://github.com/scipopt/SCIP.jl](https://github.com/scipopt/SCIP.jl)
4. Install GAMS and GAMS.jl per the instructions provided at [https://github.com/GAMS-dev/gams.jl](https://github.com/GAMS-dev/gams.jl), a community version of GAMS may be requested and is sufficient to run the included examples. In the benchmarking suite provided, BARON is accessed through the GAMS suite.
5. Run the file RSBilinear\\src\\test_benchmark.jl to generate the results. By default, the problems that were randomly generated for this trail stored in the repostory will be loaded and solved although script allows the user to randomly generated new problems if desired.

# How to Reproduce Response Surface Study Results
Follow the above instructions for configuring Julia and associated packages and run RSBilinear\\src\\rsm_model.jl.

# How to Reproduce Heat Equation Results
Follow the above instructions for configuring Julia and associated packages and run RSBilinear\\src\\heat_equation.jl.

# Associated Repositories
- [McCormick.jl](https://github.com/PSORLab/McCormick.jl) versions 0.11 or greater implements envelopes calculations for relaxations of activaiton functions.
- [EAGO.jl](https://github.com/PSORLab/EAGO.jl) versions 0.7 or greater currently includes capabilities for computing improved bilinear relaxations as detailed in the manuscript.