# EpiRegress
The code and data for the paper, EpiRegress: A method to estimate and predict the time-varying effective reproduction number, is available here. 

The data part include 4 regions: NSW for New South Wales, Australia; NZ for New Zealand; SG for Singapore; TW for Taiwan, China. 

The code consists of 5 files: 'EpiRegress_working.r' shows how to estimate the time-varying effective reproduction number (Rt) with the function 'EpiRegress' and how to plot the estimates; 'EpiRegress_rjags.r' is the rjags version of code for functions required for performing Rt estimations with EpiRegress, together with a function for plotting the Rt's; 'MCMCscript_EpiRegress.txt' is the model file for programming with rjags; 'read_X.r' loads the covariate matrix and 'read_casses.r' loads the case counts.
