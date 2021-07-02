Computer code and data for the paper "The importance of species interactions in eco-evolutionary community dynamics under climate change" by Åkesson et al., published in Nature Communications.

* Operating systems: CentOS 7, Ubuntu 19.10, Ubuntu 20.04, MacOS 10.13, Windows 10.
* Software dependencies: R (has been tested with R 4.1.0).
* Required R packages:
  - `deSolve`: Integrating differential equations.
  - `tidyverse`: Efficient data manipulation and plotting.
  - `ggpmisc`: Adding statistical information to plots.
  - `Rcpp`: Importing and compiling functions written in C/C++.
* Required non-standard hardware: none.
* Typical installation time on a normal desktop computer: no appreciable time if R is already installed. Otherwise, it is the installation time of R and the above three packages.

The following files are in the repository:

* This `README`.
* `COPYING.txt`: GNU General Public License, v3.0.
* `ecoevo_main.R`: This script runs the eco-evolutionary model. It can either be run directly from within R, or via the command prompt, by entering the following:
  `Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]`
  The input parameters are:
  - `[vbar]`: Average genetic variance (in units of degrees Celsius squared).
  - `[dbar]`: Average dispersal distance (expressed as the fraction of the pole-to-equator distance, 10,000 km).
  - `[model]`: Either `baseline`, `trophic`, `Tdep`, or `Tdep_trophic`, for the four model variants (`baseline` is without trophic interactions or temperature-dependent competition; `trophic` is with trophic interactions but no temperature-dependent competition; `Tdep` is with temperature-dependent competition but no trophic interactions; and `Tdep_trophic` is with both temperature-dependent competition and trophic interactions).
  - `[replicate]`: An integer, used to seed the random number generator with (so keeping `replicate` equal in two subsequent runs of the same model results in the exact same output).
  - `[outfile]`: Name (with path and extension) of file in which the generated data will be saved (e.g. `cc_data.csv`). The output is saved in csv format. If this argument is left blank, then the data are not saved, and instead a plot of the system's time evolution is created.
  
* `plotting_functions.R`: Functions for visualizing model results and replicating the manuscript's figures.

* `rhs_eval.cpp`: Routines for evaluating the right-hand side of the dynamical equations. They are written in C instead of R, because this is a major bottleneck in terms of speed.

* `data_10_replicates.rds`: A compressed data file containing a subsample of 10 out of our 100 replicates (to keep the file size manageable). It can be loaded in R with e.g. `dat <- readr::read_rds("data_10_replicates.rds")`. After that, one can apply the functions in `plotting_functions.R` to this dataset to reproduce the figures in our manuscript (with the caveat that the file only contains 10 out of the 100 replicates).
