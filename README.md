# Code for OGAM

### Overview

1. Directory *FNS* contains the main functions for the simulations and data applications. Specifically,
- *FNS_DataGene_Simu1* : function to generate data for simulation 1
- *FNS_DataGene_Simu2* : function to generate data for simulation 2
- *FNS_SmoBack* : functions to conduct smoothbackfitting
- *FNS_SmoBack_credit* : differs from *FNS_SmoBack* only in the stopping criterion

2. Directory *datasets* contains the datasets *credit* and *flight* for data applications. Due to the upload restrictions, we provide the preprocessed data here, which are seperated into several parts. Please **merge them before processing**. Raw data are available in the following links:
- *credit* : https://data.mendeley.com/datasets/wb3ndt69gf (P2P_Macro_Data.dta, 3 GB)
- *flight* : https://community.amstat.org/jointscsg-section/dataexpo/dataexpo2009

3. Directory *res* contains the analysis results of the simulations and data applications.

4. Directory *fig* contains the figures generated by the numerical experiments. 

5. Directory *Compare_Code* contains the code for simulation 3 to compare with other online nonparametric methods in the degenerated one-dimensional case.

### Simulations

The R scripts <font face="Candara Light">simulation1.R<font> and **simulation2.R** contain the code for the simulation study which computes the integrated mean square errors of the online and batch estimates as well as computing time, respectively. 

The R scripts **simu1_figplot.R** and **simulation2_figplot.R** contain the code to generate table XXX and Figure XXX. 

As the bandwidth selection is not our focus, all procedures of estimating the constants for bandwidths in the simulations and data applications are collected in the R script **supp-bandwidth_selection.R**. The corresponding results are stored in the directory *res*. One can load them directly when reproducing the numerical experiments.

### Data applications
#### Credit
