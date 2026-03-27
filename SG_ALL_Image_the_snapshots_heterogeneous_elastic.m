clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_heterogeneous/"));
ResultFig=fullfile('ResultFig_heterogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
sourcetype=2;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
% %The code Cp image FD result is changed for the purpose to image the
%snapshots
model=2;
figure(9);
    tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
for M=2:2:6
    SG_A_Main_code_for_imaging_the_snapshots;
end

