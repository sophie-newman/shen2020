This is a project on the bolometric quasar luminosity function.

The structure of the directory is:

codes:
	c_lib: the c code for the bolometric and extinction correction adapted from the HRH07 work
	convolution: python version of the bolometric corrections
	lf_fit: fitting codes, including the local and global fits; the output directory contains our MCMC chains
	obdata: the codes to load observational data; you might want to check the fitting code to know how to load them

data:
	all sorts of data dumped there, most of them are binned estimations of QLF
	you can find the template SED ploted in the paper, note that this SED is only for a typical luminosity (see the paper)

plots:
	codes to generate all the plots in the paper
	LFcheck: 
		various analysis on the bolometric QLF;
		you can find ways to load the QLF constrained by our global best-fit model there;
		you can also find ways to only load and plot data points there

