The whole rapid Monte Carlo simulation in the paper, including projection data simulation, scatter/attenuation correction and image reconstruction,can be achieved using this source code.

The "data" folder consists of spectrum of the incident photons, and the cross section of PMMA, Gd, Au and water. Data in this folder will be loaded at the beginning of the simualtion program.


For instance, a whole RMC simulation of Gd XFCT imaging includs following steps:
1) Run "XRF_Gd_cal.m" to obtain count of fluorescent photons;
2) Run "SCA_Gd_cal.m" to obtain count of Compton scattered photons;
3) Run "CT_Gd_sim.m" to obtain attenuation maps which will be used in system matrix calculation;
4) Run "GetSysMat_double.m" to obtain the system matrix;
5) Run "Scarate_Gd.m" to obtain the coefficients of scatter correction;
6) Run "XFCT_MLEM_Gd.m" to create a complete projcetion data and reconstruct the XFCT image using simulation data previously obtained.

Notice:
1) The simulation results in each step will be saved in "projdata" folder.
2) The parameters in each simulation program,such as the concentration and the geometry size, can be changed according to the specific requirement. Users should ensure that all the parameters in each simulation program fit the condition they want to simulate before running the program.
3) The number of particles sampled in the simulation program only affects the statiscal noise caused by the simulation program. The actual statistical noise will be simulated by adding Poisson noise after the projection data is simulated. Users should ensure that the number of particles sampled is sufficient large so that the statiscal noise caused by the simulation program is much smaller than the Poisson noise added at the end of the simulation.
4) An example of a whole RMC simulation of Gd XFCT imaging is provided.