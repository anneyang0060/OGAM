# Code for OGAM

### Overview

1. Directory *FNS* contains the main functions which are called for the simulations and data applications. Specifically,
- *FNS_DataGene_Simu1* : function to generate data for simulation setting 1.
- *FNS_DataGene_Simu2* : function to generate data for simulation setting 2.
- *FNS_SmoBack* : functions to conduct smoothbackfitting.
- *FNS_SmoBack_credit* : differs from *FNS_SmoBack* only in the stopping criterion.

2. Directory *datasets* contains the preprocessed datasets *credit* and *flight* for data applications. Raw data are available in the following links:
- *credit* : https://data.mendeley.com/datasets/wb3ndt69gf (P2P_Macro_Data.dta, 3 GB)
- *flight* : https://community.amstat.org/jointscsg-section/dataexpo/dataexpo2009

3. Directory *res* contains the analysis results of the simulations and data applications. Specifically,

- *sim1*: Numerical results of simulation setting 1, including the constants for bandwidths (batch_constants_for_bandwidths.Rdata and online_constants_for_bandwidths.Rdata), the batch estimates (batch_gam_1.Rdata--batch_gam_100.Rdata) and the online estimates  (online_gam_L'X'_1.Rdata--online_gam_L'X'_100.Rdata) with different candidate sequence lengths X=3,5,10.

- *sim2*: Numerical results of simulation setting 2, including the constants for bandwidths (batch_constants_for_bandwidths.Rdata and online_constants_for_bandwidths.Rdata), the batch estimates (batch_gam_1.Rdata--batch_gam_100.Rdata) and the online estimates  (online_gam_L'X'_1.Rdata--online_gam_L'X'_100.Rdata) with different candidate sequence lengths X=3,5,10.

- *flight*: Numerical results of the flight dataset analysis.
- *credit*: Numerical results of the credit dataset analysis.

4. Directory *fig* contains the figures generated by the numerical experiments. Specificaly,
- *theoretical_eff.pdf*: The theoretical lower bound for the relative efficiency of the proposed online estimates versus different lengths of candidate bandwidth sequences, which corresponds to Figure 4 in the main context.
- *selection_L.pdf*: The relationship between the computational cost and the length of candidate bandwidth sequence, which corresponds to Figure 5 in the main context.
- *simu_empirical_eff.pdf*: The empirical relative efficiency for the simulations that corresponds to Figure 6 in the main context.
- *simu_time.pdf*: The comparison of computing times for the simulations that corresponds to Figure 7 in the main context.
- *flights_betas.pdf*: The component function estimates for the flights dataset, which corresponds to Figure 8 in the main context.
- *credit_betas.pdf*: The component function estimates for the credit dataset, which corresponds to Figure 9 in the main context.
- *compare.pdf*: The relative efficiencies of the proposed method and other online nonparametric methods in the degenerated one-dimensional case, which corresponds to Figure 1 in the supplement.

5. Directory *Compare_Code* contains the code for simulation 3 to compare with other online nonparametric methods in the degenerated one-dimensional case.

### Workflows

#### Simulation 1 and 2

1. Make sure that the R scripts ***simulation1.R*** and ***simulation2.R*** are in the same  directory as the folder ***FNS***.
2. As the bandwidth selection is not our focus, users can choose to run Section 1 and 2 of ***supp-bandwidth_selection.R*** to generate the constants for bandwidths or to load them directly in Step 3 for convenience. If choose to load them, make sure that the following files exist:
- *res/sim1/online_constants_for_bandwidths.Rdata*
- *res/sim1/batch_constants_for_bandwidths.Rdata*
- *res/sim2/online_constants_for_bandwidths.Rdata*
- *res/sim2/batch_constants_for_bandwidths.Rdata*
3. Run the R scripts ***simulation1.R*** and ***simulation2.R*** which will load the constants for bandwidths  and call functions in the folder ***FNS*** to compute the integrated mean square errors of the online and batch estimates as well as computing times. The corresponding results will be stored in the directory ***res/sim1*** and directory ***res/sim2***. 
5. Run the R scripts ***simu_figplot.R*** to generate **Figure 4-7**. 

#### Simulation 3 in the supplementary materials

1. Run the R script ***Compare_Code/simuscript.R*** which will compute the integrated mean square errors of different online nonparametric methods in the degenerated one-dimensional case.
2. Run the R script ***Compare_Code/figplot.R*** to generate **Figure 1** of the supplement.

#### Data application 1: flight
1. Due to the upload restrictions, we provide the preprocessed data here, which are seperated into several parts. Run ***datasets/flight/flight_merge.R*** to merge them or download the raw csv data from the link in the *Overview* to the directory ***datasets/flight*** and run ***datasets/flight/flight_preprocess.R*** to preprocess the raw csv data.
2. As the bandwidth selection is not our focus, users can choose to run Section 3 of ***supp-bandwidth_selection.R*** to generate the constants for bandwidths or to load them directly in Step 3 for convenience. If choose to load them, make sure that the following files exist:
- *res/flight/online_constants_for_bandwidths.Rdata*
- *res/flight/batch_constants_for_bandwidths.Rdata*
3. Run the R script ***flight_main.R*** which will compute the online and batch estimates as well as computing times. The corresponding results will be stored in the directory ***res/flight***. 
4. Run the R script ***flight_figplot.R*** to generate **table 1** and **Figure 8** 

#### Data application 2: credit
1. Run ***datasets/credit/credit_merge.R*** to merge the seperated preprocessed data or download the raw csv data from the link in the *Overview* to the directory ***datasets/credit*** and run ***datasets/credit/credit_preprocess.R*** to preprocess the raw csv data.
2. Run Section 4 of ***supp-bandwidth_selection.R*** to generate the constants for bandwidths, or load them directly in Step 3 for convenience. If choose to load them, make sure that the following files exist:
- *res/credit/online_constants_for_bandwidths.Rdata*
- *res/credit/batch_constants_for_bandwidths.Rdata*
3. Run the R script ***credit_main.R*** which will compute the online and batch estimates as well as computing times. The corresponding results will be stored in the directory ***res/credit***. 
4. Run the R script ***credit_figplot.R*** to generate **table 2** and **Figure 9**. 
