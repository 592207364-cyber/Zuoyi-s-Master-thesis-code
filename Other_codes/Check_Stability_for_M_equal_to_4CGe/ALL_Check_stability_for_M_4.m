clearvars;clc;close all;
%% Made by Zuoyi.
sourcetype=2;
M=4;
figure(12);
tiledlayout(1,2,"TileSpacing","tight","Padding","tight");

vp0=2000.0;vs0=800.0;gamma0=vs0/vp0;
vp0=2003.0;vs0=gamma0*vp0;
  vp0=2004.0;vs0=gamma0*vp0;
  % % vp0=2005.0;vs0=gamma0*vp0;
CG_A_Main_code_for_imaging_the_snapshots;
H_Analytical_Solution_for_Stability_checking;
save('Workspace.mat');