%% This code is the main code 
%% Some codes may contain the save() at the end and load() at the beginning. 
%% No need to run all,but be careful when changing.
%% Made by Zuoyi.
FD='CG';
%% Work flow
% clearvars;close all;
% clc;
B_Set_all_the_parameters;
C_Image_the_model;
CG_D_Set_FD_coefficients;
CG_E_Check_the_stability;
if model==1
    CG_F_FD_homogeneous2D;
else
    CG_F_FD_heterogeneous2D;
end
G_Image_FD_Snapshots;

