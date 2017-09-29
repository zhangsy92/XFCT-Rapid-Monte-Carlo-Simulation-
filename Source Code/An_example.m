%-----------------------------%
%An example about a whole RMC simulation of Gd XFCT imaging including
%following steps:
%1) Run XRF_Gd_cal.m to obtain count of fluorescent photons.
%2) Run SCA_Gd_cal.m to obtain count of Compton scattered photons.
%3) Run CT_Gd_sim.m to obtain attenuation maps which will be used in 
%   system matrix calculation.
%4) Run GetSysMat_Gd_double.m to obtain the system matrix.
%5) Run Scarate.m to obtain the coefficients of scatter correction
%6) Run XFCT_MLEM_Gd.m to create a complete projcetion data and reconstruct
%   the XFCT image using simulation data previously obtained.
%-----------------------------%
run('XRF_Gd_cal.m')
run('SCA_Gd_cal.m')
run('CT_Gd_sim.m')
run('GetSysMat_Gd_double.m')
run('Scarate_Gd.m')
run('XFCT_MLEM_Gd.m')

% run('XRF_Au_cal.m')
% run('SCA_Au_cal.m')
% run('CT_Au_sim.m')
% run('GetSysMat_Au_double.m')
% run('Scarate_Au.m')
% run('XFCT_MLEM_Au.m')