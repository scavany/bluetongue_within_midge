# bluetongue_within_midge

This repository contains code to reproduce results from the following preprint:
Cavany SM, Barbera C, Carpenter M, Rodgers CP, Sherman T, Stenglein M, Mayo CM, Perkins TA **Modeling cellular co-infection and reassortment of bluetongue virus in Culicoides midges**. *bioRxiv* doi:[10.1101/2022.06.03.494504](https://doi.org/10.1101/2022.06.03.494504)

License

This code is being released under the GNU General Public License v3.0.

Overview

Contents of this repository are organized according to data, and scripts. All plots and results shown in the main text and supplementary material are generated from scripts in the scripts folder. Plots involving single infection are generated using the within_midge_barrier_based.R file. Plots involving coinfection and reassortment, including sensitivity analyses (Figs S8-9) are generated using the within_midge_barrier_based_coinfection.R file, except for those which involved varying beta (i.e. S5-7) which are generated using the within_midge_barrier_based_coinfection_beta.R file. The GAM fit to the data from Samal et al and el Hussein et al are generated using the reassortment_gam_analysis.R script. Some of the longer simulations and model fitting were run on the high-performance cluster at Notre Dame, using slightly altered versions of these scripts which are not included. The output from these runs are included as RData files in the repository to enable recreation of the manuscript's figures.
