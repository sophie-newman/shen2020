# The bolometric quasar luminosity function

---

This is a project on the bolometric quasar luminosity function. The codes are adapted from the HRH07 work by Xuejian Shen at Caltech.

- For those who are interested in the bolometric corrections and the quasar luminosity functions constrained in this work, please check out */tools/utilities.py* . You should be able to find functions that can calculate bolometric corrections and their dispersions in Section 1 of the code. 
  And find functions to return best-fit bolometric QLFs or predicted QLFs in bands in Section 2. Remember to compile the c code in */tools/clib/* before use.

- For those who want to use the observational data used in this work, please check the fitting code */codes/lf_fit/lf_fitter.py* or the plot codes in */plots/LFcheck/plot_allbands.py*. There, you should be able to find a function that loads all the observational data in */codes/obdata/*.

# For those who are interested in more details, the structure of the directory is as follows:

---

## codes:
* bolometric_correction: 
	* codes to calculate the bolometric correction based our SED model; 
	* the ensemble of SEDs to measure the dispersion of the bolometric corrections are also saved there
* c_lib: the c code for the bolometric and extinction correction adapted from the HRH07 work
* convolution: python version of the bolometric corrections
* lf_fit: 
	* fitting codes, including the local and global fits; 
	* the output directory contains our MCMC chains
* obdata: 
	* the codes to load the observational data in individual paper; 
	* you might want to check the fitting code to know how to load the data from all the papers

---

## data:
all sorts of data dumped there, most of them are binned estimations of the QLF (check the fitting code or the plot codes to know how to load them)

you can find the template SED ploted in the paper (MySED.dat), note that this SED is only for a typical luminosity (see the paper)

---

## plots:
codes to generate all the plots in the paper

some important plots are:

* LFcheck: 
	* various analysis on the bolometric QLF;
	* you can find ways to load the QLF constrained by our global best-fit model there;
	* you can also find ways to only load and plot data points there
* ionization: predictions on the ionization budget contributed by quasars and predictions on the Lyman-continuum emissivity of quasars
* integrated background: CXB predictions
* Fit_parameters: analysis on the evolution of the DPL LF parameters

---

If you have any questions, please contact: xshen@caltech.edu

