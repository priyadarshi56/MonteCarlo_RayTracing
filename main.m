%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   MAIN file for MC ray-tracing simulation in two dimensional domain
%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; % close all

DateTimeStart=datestr(now,'dd-mm-yyyy_HH-MM-SS');
constants;
input_parameters;
input_geometry;
meshing_domain;
Poisson_solver;
%----- Multiplication factor ----------------------------------------------
tic
[MF_avg]=MF_calculator(Lx,dLx,Ly,dLy);
MF_time=toc;
fprintf('\nMF time= %g seconds.\n',MF_time);
%----- MC Simulation ------------------------------------------------------
run_raytracing;
% run_raytracing_parallel
%----- TE Parameter calculation and plotting ------------------------------
calculate_TE_parameters;
% plot_domain_meshgrid
% plot_flux_TDF;
% plot_TE_parameters
%--------------------------------------------------------------------------
Sim_time=mesh_time+MF_time+poisson_time+MC_time_taken;
fprintf('\nTotal simulation time is %g seconds.',Sim_time);
fprintf('\n------------------------- Simulation Ends -------------------------\n')
DateTimeEnd=datestr(now,'dd-mm-yyyy_HH-MM-SS');
save_files;

write;




