# Repository for the paper "Within-host viral growth and immune response rates predict FMDV transmission dynamics for African Buffalo"

Folders in this repository are: <br />
  1. Archive
  2. Finalized fitting
  3. Output analysis 

The <ins> Archive </ins> folder contains previous versions of the model fitting, including efforts to fit time series data, such as SAA data for innate immune response, which was not part of the finalized analysis.  Most users will want to look instead at the <ins> Finalized fitting </ins> and <ins> Output analysis </ins> folders.

## Finalized fitting
The <ins> Finalized fitting </ins> folder contains the (i) finalized code for model fitting as described in the methods section and SI of the paper, (ii) the input experimental data formatted to work with this code, (iii) the code to generate parameter confidence intervals for both the needle and contact-infected hosts, (iv) the code to generate the time plots shown in figures 2 and S3 of the manuscript, (v) output parameters in .mat format.  Specific files in this folder (with brief description) are:

### Data files
- <ins>HaptoDataContact.csv</ins> contains the Haptoglobin time series data (a measure of innate immune response) for each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are micrograms per mL (log base 10 is taken prior to fitting).
- <ins>VNTContactData.csv</ins> contains the Viral neutralization (VNT) time series data for each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are log base 10 VNT.
- <ins>ViremiaContact.csv</ins> contains the Viral load time series data afor each of the twelve contact-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since contact was initiated).  Units are log base 10 genome copies per mL.
- <ins>HaptoDataNeedle.csv</ins> contains the Haptoglobin time series data (a measure of innate immune response) for each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are micrograms per mL (log base 10 is taken prior to fitting).
- <ins>VNTNeedleData.csv</ins> contains the Viral neutralization (VNT) time series data for each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are log base 10 VNT.
- <ins>ViremiaNeedle.csv</ins> contains the Viral load time series data afor each of the twelve needle-infected hosts.  Rows correspond to each host, and columns to observation day (in terms of time since experiment start).  Units are log base 10 genome copies per mL.### Program files
- <ins>InfectionStartTImes.csv</ins> contains the posterior infection start time estimates generated in our previous work [LINK]
- <ins>params_noise.mat</ins> is a matlab matrix file containing the 10,000 paramter estimates per host used to generate figure 2 of the manuscript
- <ins>params_noise_needle.mat</ins> is the matlab matrix file containing the 10,000 paramter estimates per host used to generate figure S3 of the manuscript

## Program files
- <ins>ContactDataHaptoStimes.m</ins> is the matlab code to generate a single instance of the Monte-Carlo simulations for a specified contact-infected host
- <ins>NeedleDataHaptoStimes.m</ins? is the matlab code to generate a single instance of the Monte-Carlo simulations for a specified needle-infected host
## Output analysis

### Data files

### Program files
