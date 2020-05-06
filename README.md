Computer code for the manuscript "The importance of species interactions in spatially explicit eco-evolutionary community dynamics under climate change".

* Operating systems: CentOS 7, Ubuntu 19.10, MacOS 10.13, Windows 10
* Software dependencies: R 3.6.1 (software has been tested on R 3.6.1 and above)
* Required R packages:
  - `deSolve` (for integrating differential equations)
  - `tidyverse` (for efficient data manipulation and plotting)
* Required non-standard hardware: none
* Typical install time on a normal desktop computer: no appreciable time if R is already installed. Otherwise, it is the installation time of R.

There are two R scripts in this repository:

* `ecoevo_main.R`: This script runs the eco-evolutionary model. It can either be run directly from within R, or via the command prompt, by entering the following:
  `Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]`
  The input parameters are:
  - `[vbar]`: average genetic variance (in units of degrees Celsius squared)
  - `[dbar]`: average dispersal distance (expressed as the fraction of the pole-to-equator distance, 10,000 km)
  - `[model]`: either `baseline`, `trophic`, `Tdep`, or `Tdep_trophic`, for the four model variants (`baseline` is without trophic interactions or temperature-dependent competition; `trophic` is with trophic interactions but no temperature-dependent competition; `Tdep` is with temperature-dependent competition but no trophic interactions; and `Tdep_trophic` is with both temperature-dependent competition and trophic interactions)
  - `[replicate]`: an integer, used to seed the random number generator with (so keeping `replicate` equal in two subsequent runs of the same model results in the exact same output)
  - `[outfile]`: name (with path and extension) of file in which the generated data will be saved (e.g. `cc_data.csv`). The output is saved in csv format. If this argument is left blank, then the data are not saved, and instead a plot of the dynamics is created.
  
* `plotting_functions.R`: functions for visualizing model results and replicating the figures of the manuscript

