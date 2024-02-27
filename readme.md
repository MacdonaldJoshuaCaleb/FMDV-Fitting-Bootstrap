# Repository for the paper "Within-host viral growth and immune response rates predict FMDV transmission dynamics for African Buffalo"

## Abstract
Infectious disease dynamics operate across biological scales: pathogens replicate within hosts but transmit among populations. Functional changes in the pathogen-host interaction thus generate cascading effects across organizational scales. We investigated within-host dynamics and among-host transmission of three strains (SAT1, 2, 3) of foot-and-mouth disease viruses (FMDVs) in their wildlife host, African buffalo. We combined data on viral dynamics and host immune responses with mathematical models to ask (i) How do viral and immune dynamics vary among strains?; (ii) Which viral and immune parameters determine viral fitness within hosts?; and (iii) How do within-host dynamics relate to virus transmission? Our data reveal contrasting within-host dynamics among viral strains, with SAT2 eliciting more rapid and effective immune responses than SAT1 and SAT3. Within-host viral fitness was overwhelmingly determined by variation among hosts in immune response activation rates but not by variation among individual hosts in viral growth rate. Our analyses investigating across-scale linkages indicate that viral replication rate in the host correlates with transmission rates among buffalo and that adaptive immune activation rate determines the infectious period. These parameters define the virus's relative basic reproductive number ($\mathcal{R}_0$), suggesting that viral invasion potential may be predictable from within-host dynamics.

## Authors
Joshua C. Macdonald, Hayriye Gulbudak, Brianna Beechler, Erin E. Gorsich, Simon Gubbins, Eva Perez-Martin, and Anna E. Jolles

## Corresponding authors contact
- AJ: aejolles-AT-gmail.com, HG: h.gulbudak-AT-louisiana.edu
  
## Contact for questions about program files
- JCM: joshuamac-AT-tauex.tau.ac.il (temporary, postdoc email), or macdonald.j.caleb-AT-gmail.com (stable)
  
## The files
Folders in this repository are: <br />
  1. Finalized fitting
  2. Output analysis 

## Finalized fitting
The <ins> Finalized fitting </ins> folder contains the (i) finalized code for model fitting as described in the methods section and SI of the paper, (ii) the input experimental data formatted to work with this code, (iii) the code to generate parameter confidence intervals for both the needle and contact-infected hosts, (iv) the code to generate the time plots shown in figures 2 and S3 of the manuscript, (v) output parameters in 3D (host x paramter x simulation number) .mat format. needed to run the CI generation and plot generation files (vi) For useers that may not have access to matlab there are 2D csv versions of these files included.  Specific files in this folder (with brief description) are:

### Data files
- <ins>HaptoDataContact.csv</ins> contains the Haptoglobin time series data (a measure of innate immune response) for each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are micrograms per mL (log base 10 is taken prior to fitting).
- <ins>VNTContactData.csv</ins> contains the Viral neutralization (VNT) time series data for each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are log base 10 VNT.
- <ins>ViremiaContact.csv</ins> contains the Viral load time series data afor each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are log base 10 genome copies per mL.
- <ins>HaptoDataNeedle.csv</ins> contains the Haptoglobin time series data (a measure of innate immune response) for each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are micrograms per mL (log base 10 is taken prior to fitting).
- <ins>VNTNeedleData.csv</ins> contains the Viral neutralization (VNT) time series data for each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are log base 10 VNT.
- <ins>ViremiaNeedle.csv</ins> contains the Viral load time series data afor each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are log base 10 genome copies per mL.### Program files
- <ins>InfectionStartTImes.csv</ins> contains the posterior infection start time estimates generated in our previous work (https://doi.org/10.5281/zenodo.5121203) 
- <ins>params_noise.mat</ins> is the 3D matlab matrix file containing the 10,000 paramter estimates per host used to generate figure 2 of the manuscript
- <ins>params_noise.csv</ins> is the 2D csv version of the matlab matrix file of the same nae to allow access to indivdual host parameter estimates for those without matlab
- <ins>params_noise_needle.mat</ins> is the 3D matlab matrix file containing the 10,000 paramter estimates per host used to generate figure S3 of the manuscript
- <ins>params_noise_needle.csv</ins> is the 2D csv version of the matlab matrix file of the same nae to allow access to indivdual host parameter estimates for those without matlab
- <ins> QualCompInd.csv</ins> has median paramter estimates for both within and between-host parameters across all contact infected hosts
  
### Program files
#### Support files that don't require end user modification 
- <ins>ContactDataHaptoStimes.m</ins> is the matlab code to generate a single instance of the Monte-Carlo simulations for a specified contact-infected host
- <ins>NeedleDataHaptoStimes.m</ins> is the matlab code to generate a single instance of the Monte-Carlo simulations for a specified needle-infected host
- <ins>GetData.m</ins> is the file to generate a new time series given a model fit
- <ins>ViralGrowth.m</ins> generates initial estimates of the initial viral load and viral growth rate as described in the manuscript
- <ins>GetDataInitial.m</ins> generates new inital data points for immune response data given a drawn infection start time as described in the manuscript
#### Scripts with function calls to replicate figures and confidence interval generation 
- <ins> tplotsNew.m</ins> generates the time plots presented in figure 2
- <ins> tplotsNewNeedle.m</ins> generates the time plots presented in figure S3

#### Script go generate 10,000 parameter estimates per host 
- <ins>GetParams.m</ins> generates the paramter estimates over 10,000 simulations for all hosts (both needle and contact infected) by calling support function files 

## Output analysis
<ins> Output analysis</ins> contains the mean parameter estimates for each serotype over the 10,000 Monte-Carlo simulations per host as well as the python script that generates the correlation and regression analyses as well as figures 3-5 of the manuscript.

### Data files
- <ins>FMDVMeansMay9.csv</ins> 10,000 sample means per serotype, columns are individual parameters, rows are individual fitting interations
- <ins> QualCompInd.csv</ins> has median paramter estimates for both within and between-host parameters across all contact infected hosts
### Program files
- <ins>Fitting_analysis.ipynb</ins> Python notebook for correlation and regression analyses as well as generation of figures 3-5.
