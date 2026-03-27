clearvars;clc;close all;
%% Made by Zuoyi.
sourcetype=1;
M=6;
figure(12);
tiledlayout(1,2,"TileSpacing","tight","Padding","tight");

vp0=2876.2;
vp0=2876.8;
CG_A_Main_code_for_imaging_the_snapshots;
H_Analytical_Solution_for_Stability_checking;
save('Workspace.mat');