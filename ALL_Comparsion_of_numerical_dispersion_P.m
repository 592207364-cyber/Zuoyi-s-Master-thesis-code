clearvars;clc;close all;
addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_homogeneous/"));
ResultFig=fullfile('ResultFig_homogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));
color = [
    0.620,  0.792,  0.882;   % CGFD Elastic 1 (Light Blue)
    0.420,  0.682,  0.839;   % CGFD Elastic 2 (Medium Light Blue)
    0.192,  0.510,  0.741;   % CGFD Elastic 3 (Medium Dark Blue)
    0.031,  0.318,  0.612;   % CGFD Elastic 4 (Dark Blue)
    0.843,  0.098,  0.110;   % CGFD Acoustic  (Brick Red)
    0.000,  0.000,  0.000;   % SGFD           (Pure Black)
    0.500,  0.500,  0.500;   % Padding (Gray)
    0.500,  0.500,  0.500;   % Padding (Gray)
    0.500,  0.500,  0.500    % Padding (Gray)
];
for phi=0:pi/4:pi/4
    figure(4+phi*4/pi);
    tiledlayout(3,1,"TileSpacing","tight","Padding","tight");
    for M=2:2:6
        figure(4+phi*4/pi);
        nexttile(M/2);hold on;
        
        model=1;
         
        
        sourcetype=2;%=2: Vertical force;
        colorindex=1;
        for gamma0=0.2:0.2:0.8
            CG_A_Main_code_for_numerical_dispersion;
            plot(kh(1:314-2*(phi==0)-gamma0*10*(phi==0)),cgcfdvp_c(1:314-2*(phi==0)-gamma0*10*(phi==0)),Color=color(colorindex,:),LineWidth=2,DisplayName=['Elastic CG-FD P-wave dispersion, γ=',num2str(gamma0)]);hold on;grid on;
            colorindex=colorindex+1;
        end
        sourcetype=1;% =1: Pressure;
         CG_A_Main_code_for_numerical_dispersion;
         plot(kh(1:314-2*(phi==0)-10*(phi==0)),cgcfdvp_cacoustic(1:314-2*(phi==0)-10*(phi==0)),Color=color(5,:),LineWidth=2,DisplayName='Acoustic CG-FD P-wave dispersion',Visible='on');hold on;
        
         sourcetype=2;%=2: Vertical force;
         SG_A_Main_code_for_numerical_dispersion;
          plot(kh(1:314-2*(phi==0)-12*(phi==0)*(M==2)-2*(phi==pi/4)*(M==2)),sgcfdvp_c(1:314-2*(phi==0)-12*(phi==0)*(M==2)-2*(phi==pi/4)*(M==2)),'k',LineWidth=2,DisplayName='SG-FD P-wave dispersion',Visible='on');hold on;grid on;
        lgd = legend('Location', 'southwest');
        lgd.ItemTokenSize = [15, 18];
        legend show;
        % legend('SG-v','CG-acoustic','CG-u_x γ=0.2','CG-u_x γ=0.4','CG-u_x γ=0.6','CG-u_x γ=0.8',location='southwest');hold on;grid on;
        % legend(FontSize=20);
        xlim([0 pi]);
        ylim([0.62 1.03]);xlabel('kh');ylabel('c_f_d/c');
        title_str = sprintf('c_{fd}/c - kh, \\theta = %.0f^\\circ, 2M = %d, v_p = %.1f m/s', ...
                    phi/pi*180, M, vp0);
        title(title_str, 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
        set(gcf, 'Position', [0, 0, 600, 900], 'Color', 'w');
    end
    saveas(gcf,fullfile([ResultFig,'/Numerical_Dispersion'],['cfd_c_kh_theta_',num2str(phi/pi*180),'_2M_2_4_6.png']));
end