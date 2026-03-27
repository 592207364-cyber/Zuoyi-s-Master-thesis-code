clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_homogeneous/"));
ResultFig=fullfile('ResultFig_homogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
sourcetype=2;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
model=1;
% %The code Cp image FD result is changed for the purpose to image the
%snapshots
% figure(9);
% tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
for M=2:2:6
    CG_A_Main_code_for_imaging_the_snapshots;
    save('Workspace.mat');
    %I_Show_the_anime;
end

