# Balancing Inferential Integrity and Disclosure Risk Via Model Targeted Masking and Multiple Imputation
# Author Contributions Checklist Form

## Data

**Abstract (Mandatory)**

The data consist of scleroderma patient registry data collected by the Canadian Scleroderma Research Group (CSRG) under the direction of Dr. Murray Baron from McGill University. Patients in the registry are recruited from 15 centers across Canada. Our analyses included available patient data from the baseline and first follow-up visits of patients with baseline and first follow-up visits between September 2004 and May 2016.

**Availability (Mandatory)**

Due to privacy concerns and the nature of the patient consent agreements, the original data cannot be made available to the public.  Researchers who want to obtain access to the data should contact one of the authors of the paper (Dr. Steele) for details on the procedures for requesting access to the data.

**Description (Mandatory if data available)**

For reproducibility, we have created a pseudo data set (&quot;csrg.csv&quot;) that mimics the real CSRG data set used in the paper. This pseudo CSRG data set is available in the on-line supplementary material.

**Optional Information (complete as necessary)**

N/A

## Code

**Abstract (Mandatory)**

The C++ and R code to reproduce the main results using the pseudo data set (Tables 2-4 in the main paper; the corresponding results using the pseudo data set are reported in Web Tables 17-19 in the online supplementary material) and the main simulation results provided in the online supplementary material. Specifically, the file simu.zip contains code to re-run the simulation study and reproduce the main simulation results reported in Web Figures 2-8 in the online supplementary material. The file real.zip contains the pseudo CSRG data set and the code to reproduce the data analysis results using this pseudo data set reported in Web Tables 17-19 in the supplementary material.

**Description (Mandatory)**

The real.zip file contains all necessary code to reproduce the main results (Web Tables 17-19) using the pseudo CSRG data set:

- The pseudo CSRG data set: csrg.csv
- joint\_multivariate\_sensitive\_as\_covariate.cpp contains C++ functions to produce synthetic data sets using the proposed MI-DA method.
- run\_synthpop.r contains R code to impute the missing data (by calling csrg\_MI.r with the seed set); generate synthetic data sets for the MI complete data sets using the CART and Norm+Logit approaches implemented in the synthpop package (by calling syn\_synthpop.r) and produce disclosure risks and data utilities in both the work disability and ILD analysis (by calling utility functions defined in combine.r).
- syn\_synthpop.r contains the R code to generate synthetic data sets using the CART and Norm+Logit approaches implemented in the synthpop package.
- run\_DA-MI.r contains R code to impute the missing data (by calling csrg\_MI.r with the seed set); generate synthetic data sets for the MI complete data sets using the the proposed DA-MI approach (by calling generate\_w.r to create copies of the pseudo variables in the data augmentation step and calling syn\_DA-MI.r) and produce disclosure risks and data utilities in both the work disability and ILD analysis (by calling defined utility functions defined in combine.r).
- syn\_DA-MI.r contains R code to generate synthetic data sets using the proposed DA-MI approach, which involves calling the joint\_multivariate\_sensitive\_as\_covariate.cpp.
- combine.r contains R functions to impute missing data, obtain the data utility measures and disclosure risks for the synthetic data sets.
- disclosure\_risk.r contains R code to calculate different data utility measures and disclosure risks using two strategies for the synthetic data sets, where the risk calculation is adapted based on the R code found at

https://www2.stat.duke.edu/~jerry/JPSMCourseR/programs/calculate\_risks.R.

- summarize.r contains R code to generate the &quot;.tex&quot; files for Web Tables 17-19.

The simu.zip file contains all necessary code to reproduce the main simulation results (Web Figures 2-8):

- joint\_multivariate\_sensitive\_as\_covariate.cpp contains C++ functions to produce synthetic data sets using the proposed MI-DA method.
- gen\_data.r contains R functions to generate simulated data sets.
- generate\_synthetic\_synthpop.r contains R functions to generate simulated datasets by calling gen\_data.r and setting the seed as the data set number, produce synthetic data sets using Norm+Logit and CART methods implemented in synthpop package in R, and obtain the data utility measures and disclosure risks for the synthetic datasets, by calling the functions defined in fun\_for\_summary\_simu.r.
- generate\_synthetic.r contains R code to generate simulated datasets by calling gen\_data.r and setting the seed as the data set number, and provides R interface to call joint\_multivariate\_sensitive\_as\_covariate.cpp to produce synthetic data sets using the proposed MI-DA method.
- fun\_for\_summary\_simu.r contains supporting R functions to summarize the simulation results.
- for\_summary.r contains specific R functions to calculate different data utility measures and disclosure risks using two strategies for the synthetic data sets, where the risk calculation is adapted based on the R code found at

https://www2.stat.duke.edu/~jerry/JPSMCourseR/programs/calculate\_risks.R.

- summary\_DA.r contains R functions to obtain the data utility measures and disclosure risks for the synthetic data sets obtained from the proposed MI-DA method, by calling the functions defined in fun\_for\_summary\_simu.r and fun\_for\_summary\_simu.r.
- function\_for\_plot.r contains supporting R functions to generate figures.
- figures.r contains the functions to create the Web Figures 2-8 by calling function\_for\_plot.r.

**Optional Information (complete as necessary)**

The C++ code was developed using Open Source C++ library Scythe-1.0.3, please refer to [http://scythe.lsa.umich.edu](http://scythe.lsa.umich.edu/) for more details and installation instructions.

The version of the R software used is R 3.6.0.

The following R packages are used: synthpop (version 1.5-1), mi (version 1.0), Hmisc (version 4.3-0).

## Instructions for Use

**Reproducibility (Mandatory)**

To reproduce the main data analysis results (Web Tables 17-19) in the online supplementary material, run the R files in the following order:

1. Run run\_DA-MI.r to load the data file csrg.csv and obtain the final result files by using the proposed DA-MI approach:

result\_csrg\_Age\_std=1\_Gender=6\_NewWhite=6\_MoreThanHighSchool=6\_K=20\_r=3\_3.Rdata

and

result\_csrg\_Age\_std=1\_Gender=6\_NewWhite=6\_MoreThanHighSchool=6\_K=20\_r=3\_5.Rdata

using the final tuning setting considered for the real CSRG data set.

Note, calling run\_DA-MI.r will also generate a benchmark result file csrg\_benchmark.Rdata used for comparison in Web Tables 17-19.

1. Run run\_synthpop.r to load the same csrg.csv data file and obtain the final result files by using the Norm+Logit and CART approaches:

synthpop\_d=3.Rdata and synthpop\_d=5.Rdata

Note, calling run\_synthpop.r generates the same benchmark result file csrg\_benchmark.Rdata.

1. Run summary.r to load the above result files (obtained from Steps 1 and 2) and generate the following latex tables:

csrg\_simu\_risk.tex for Web Table 17

csrg\_simu\_est.tex for Web Table 18

csrg\_simu\_sd\_overlap.tex for Web Table 19.

To reproduce the main simulation results (Web Figures 2-8) in the supplementary material, run the R files in the following order:

1. Run generate\_synthetic.r to generate simulated data sets and create multiple synthetic data sets using our proposed MI-DA method under different settings of the factors (i.e., tuning parameters and proportion of synthesis).
2. Run generate\_synthetic\_synthpop.r to generate the same simulated data sets as in 1, create multiple synthetic data sets using Norm+Logit and CART methods, and evaluate their data utility measures and disclosure risks.
3. Run summary\_DA.r to obtain the data utility measures and disclosure risks in the synthetic data sets obtained using the proposed MI-DA method.
4. Run figures.r to generate Web Figures 2-8.

**Replication (Optional)**

N/A
