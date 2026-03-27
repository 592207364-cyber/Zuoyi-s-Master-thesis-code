clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_homogeneous/"));
ResultFig=fullfile('ResultFig_homogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
model=1;
sourcetype=1;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
FD='CG';
figure(6);
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

% pause;
clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_homogeneous/"));
ResultFig=fullfile('ResultFig_homogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
model=1;
sourcetype=1;% =1: Pressure; =2: Vertical force;  =3: Horizontal force
FD='SG';
figure(6);
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


