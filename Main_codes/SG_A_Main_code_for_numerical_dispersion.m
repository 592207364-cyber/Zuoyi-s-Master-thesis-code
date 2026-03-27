%% This code is the main code 
%% Some codes may contain the save() at the end and load() at the beginning. 
%% No need to run all,but be careful when changing.
%% The code is based on the code of Thomas.Bohlen but sorted and 
%% changed by Zuoyi to make the better comparsion.
FD='SG';
%% Work flow
% clearvars;close all;
%clc;
B_Set_all_the_parameters;
SG_D_Set_FD_coefficients;
SG_E_Check_the_stability;
SG_F_FD_2D;
SG_E_homogeneous_cfd_c;
