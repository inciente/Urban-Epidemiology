README File for software and data supplementary to "Understanding the role of urban design in disease spreading"

This document describes the order in which files are run, as well as what comes out of each one.

1. Run /Demography/onto_grid.m 
	This file loads census data from 'jalisco_2010_ageb_data.csv', extracts variables relevant to our study and then interpolates their information onto the triangular grid defined in 'malla892nod.mat'.
	At the end, this assigns population density values for children and adults in each grid element in the file 'city.mat'.
	
2.  Run /Demography/get_TAS.m
	
	This uses data from the DENUE economic census (found in 'denue_inegi_2016_redux.csv') to calculate TAS.
	Every type of establishment is associated with a code, which are explained in 'TAS_activity_codes.csv', then weighted using coefficients TAS_coef.
	This way, TAS is calculated for every establishment registered in the census, and the total TAS within each grid element is then summed. 
	We also use the activity code for health care facilities to get a database of their locations. 
	
	Later on, we load SNIE data on the location and enrollment of all educational institutions in the city. These are the files 'prepasGDL.csv','unisGDL.csv', and 'escuelasGDL.csv'.
	Once we have this, we repeat the procedure of weighting TAS (by enrollment) within each grid element and we save that within city.mat. 

3. Calibrate and run simulations in setup_simulation.m

	Make sure that you have  GPU setup with CUDA to run this, or else erase the gpuArray commands and run in your CPU (slower).
	This is where the simulations are run. Begin by setting all mobility and disease parameters, select the number of initial conditions you want to run and how long your simulations should be. 
	TAS data will be used to generate origin destination matrices and produce the 'mobility.mat' file (invokes get_OD_matrix function).
	
	
Extra: want to modify the structure and agglomeration of activities in GDL?
	Run sensibility.m 
		This script takes in a series of exponentials (<1 reduces Gini coefficients, > increases them) to transform the spatial distribution of housing and daily activities in GDL.
		Look into the script "alternate_geography.m" to see exactly how this is done. You can choose whether to preserve the spatial structure of GDL or scramble things up.
		In the end you will have Gini coefficients, R_0 and epidemic size for the hypothetical cities you ran simulations on. 

Don't hesistate in reaching out to Noel G. Brizuela (nogutier@ucsd.edu) if you have any questions. 

