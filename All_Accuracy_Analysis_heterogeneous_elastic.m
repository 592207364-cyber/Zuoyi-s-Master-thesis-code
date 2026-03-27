clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_heterogeneous/"));
ResultFig=fullfile('ResultFig_heterogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
model=2;
sourcetype=2;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
FD='CG';
figure(6);
tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
figure(7);
tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
% for M=2:2:6
for M=2:2:6
    disp(['-------------CG,2M=',num2str(M)])
    B_Set_all_the_parameters;
   
    CG_D_Set_FD_coefficients;
    CG_E_Check_the_stability;
    CG_F_FD_homogeneous2D;
    H_Analytical_Solution;

end
save('Workspace.mat');
% pause;
clearvars;clc;close all;
tempData1 = load('Workspace.mat', 'routputx'); 
tempData2 = load('Workspace.mat', 'routputz'); 
routputx1 = tempData1.routputx;
routputz1 = tempData2.routputz;
clear tempData1;
clear tempData2;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_heterogeneous/"));
ResultFig=fullfile('ResultFig_heterogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
model=2;
sourcetype=2;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
FD='SG';
figure(6);
tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
figure(7);
tiledlayout(1,3,"TileSpacing","tight","Padding","tight");
% for M=2:2:6
for M=2:2:6
    disp(['-------------SG,2M=',num2str(M)])
    B_Set_all_the_parameters;
    SG_D_Set_FD_coefficients;
    SG_E_Check_the_stability;
    SG_F_FD_2D;
    H_Analytical_Solution;
end


