clearvars;clc;close all;

CG_EE_Stability_Grid_Scan;

addpath(genpath("Main_codes/"));
addpath(genpath("ResultFig_homogeneous/"));
ResultFig=fullfile('ResultFig_homogeneous');
% rmpath(genpath("Main_codes/"));
% rmpath(genpath("ResultFig/"));

MM=2:2:20;
rmaxcg_p=zeros(1,length(MM));
 rmaxcg_pp=zeros(1,length(MM));
rmaxsg=zeros(1,length(MM));
for j=1:length(MM)
    M=MM(j);
    CG_D_Set_FD_coefficients;
    SG_D_Set_FD_coefficients;
    sumabsbm=0;
    sumbm=0;
    sumam=0;
    sumabsam=0;
    sumabsb_m=0;
    for m=1:M/2
        sumam=sumam+am(m);
        sumabsam=sumabsam+abs(am(m));
        sumabsb_m=sumabsb_m+abs(b_m(m));
        sumabsbm=sumabsbm+abs(bm(m));
        sumbm=sumbm+bm(m);
    end
    rmaxcg_p(j)=sqrt(2)/sqrt(2*sumabsam+2*sumam);
    
    rmaxsg(j)=1/sqrt(2)/sumabsb_m;
end
figure(8);
%% Don't Delete this plot(this plot is for PPT presentation)%%Cubic SPline Interpolation
%rmaxcg_p=RMAXCG-rmaxcg_p;%this line code is just to calculate the difference between gird scan result and von neumann result.
xi=linspace(min(MM),max(MM),300);yi=interp1(MM,rmaxcg_p,xi,'spline');yii=interp1(MM,rmaxsg,xi,'spline');
  % plot(xi,yi,'HandleVisibility','off',LineWidth=1,Color='k');hold on;grid on;
%  plot(xi,yii,'HandleVisibility','off',LineWidth=2,Color='k');hold on;grid on;
%%





 % plot(MM,rmaxcg_p,'o-',DisplayName='CG-acoustic von Neumann analysis result',Color=color(5,:));hold on;grid on;
  plot(MM,rmaxsg,'ko-',DisplayName='SG-FD');hold on;
title('Stability limits of 2D CG-FD and SG-FD')
%title('Stability limits of 2D CG-FD-Elastic')
% title({'Stability limits of 2D Acoustic CG-FD:', 'Grid Scan method vs. von Neumann analysis'}, 'FontSize', 12, 'Interpreter', 'none')
xlabel('Accuracy order 2M')
ylabel('r_m_a_x');
axis tight;
% plot([4,4],[0.8012,0.80152],'r*');
legend show;

%% box
% axes('Position', [0.35, 0.55, 0.2, 0.2]); 
% box on; 
% plot(x,y,'o-','DisplayName','CG-acoustic grid scan result',color=color(5,:));hold on;grid on;
% % plot(xi,yi,'HandleVisibility','off','LineWidth',1,Color=color(5,:));hold on;grid on;
% hold on; grid off;
% % plot(xi,yi,'HandleVisibility','off',LineWidth=1,Color='k');hold on;grid on;
% plot(MM,rmaxcg_p,'ko-',DisplayName='CG-acoustic von Neumann analysis result');hold on;grid on;
% xlim([5.98, 6.02]);xl = xlim;xticks( ceil(xl(1)/2)*2 : 2 : floor(xl(2)) );
% ylim([0.5750, 0.5754]);
% set(gca, 'Color', 'w', 'FontSize', 10, 'FontName', 'Times New Roman');
% title('Zoomed-in Detail1', 'FontSize', 10);
% 
% axes('Position', [0.65, 0.55, 0.2, 0.2]); 
% box on; 
% plot(x,y,'o-','DisplayName','CG-acoustic grid scan result',color=color(5,:));hold on;grid on;
% % plot(xi,yi,'HandleVisibility','off','LineWidth',1,Color=color(5,:));hold on;grid on;
% hold on; grid off;
%  % plot(xi,yi,'HandleVisibility','off',LineWidth=1,Color='k');hold on;grid on;
% plot(MM,rmaxcg_p,'ko-',DisplayName='CG-acoustic von Neumann analysis result');hold on;grid on;
% xlim([7.98, 8.02]);xl = xlim;xticks( ceil(xl(1)/2)*2 : 2 : floor(xl(2)) );
% ylim([0.5546, 0.5550]);
% set(gca, 'Color', 'w', 'FontSize', 10, 'FontName', 'Times New Roman');
% title('Zoomed-in Detail2', 'FontSize', 10);
%% Paper
% saveas(gcf,fullfile([ResultFig,'/Rmaxlimit'],'rmax_2M.png'));
%% PPT
saveas(gcf,fullfile([ResultFig,'/Rmaxlimit'],'rmax_2M_for_PPT.png'));
% plot(zeros(1,length(gammarange))+2,1./(sqrt(1+gammarange.^2)),'r+',LineWidth=2,DisplayName='1/sqrt(1+γ^2)');hold on;grid on;
saveas(gcf,fullfile([ResultFig,'/Rmaxlimit'],'rmax_2M_for_PPT_with_checking.png'));