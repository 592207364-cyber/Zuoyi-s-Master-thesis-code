%% Analytical solution
%% GF2Del   % % [vx,vy,t]=GF2Del(fc,vp0,vs0,rho,dt,T,x,y)

switch FD
    case 'CG'
        char1='u';
        as_ux=zeros(nr,nt);as_uz=zeros(nr,nt);%% as means analytical solution
        as_p=zeros(nr,nt);
        R=sqrt((xr-xs).^2+(zr-zs).^2);
        for c=1:nr
            [as_ux(c,:),as_uz(c,:),~]=Z_Function_GF2Del(fc,vp0,vs0,rho0,dt,T,(xr(c)-xs)*dh,(zr(c)-zs)*dh,FD);
            %  ux(c,:)=ux(c,:)+c;
            % uz(c,:)=uz(c,:)+c;
            [as_p(c,:),t]=Z_Function_GF2Dac(fc,vp0,R(c),dt,T);
        end

        for c=1:nr
            routputx(c,:)=ux(xr(c),zr(c),:);
            routputz(c,:)=uz(xr(c),zr(c),:);
            routputp(c,:)=p(xr(c),zr(c),:);
        end
    case 'SG'
        char1='v';
        %% Analytical solution
        as_vx=zeros(nr,nt);as_vz=zeros(nr,nt);%% as means analytical solution
        as_p=zeros(nr,nt);
        %give the real location of pressure,vx,vz
        xr_p=xr*dh;zr_p=zr*dh;xs_p=xs*dh;zs_p=zs*dh;
        xr_vx=xr*dh+0.5*dh;zr_vx=zr*dh;
        xr_vz=xr*dh;zr_vz=zr*dh+0.5*dh;
        %calculate the real distance of receivers and source of pressure points
        R=sqrt((xr_p-xs_p).^2+(zr_p-zs_p).^2);
        for c=1:nr
            [as_vx(c,:),as_vz(c,:),~]=Z_Function_GF2Del(fc,vp0,vs0,rho0,dt,T,(xr(c)-xs)*dh,(zr(c)-zs)*dh,FD);
            %  vx(c,:)=vx(c,:)+c;
            % vz(c,:)=vz(c,:)+c;
            [as_p(c,:),t]=Z_Function_GF2Dac(fc,vp0,R(c),dt,T);
        end
        as_ux=as_vx;as_uz=as_vz;
        routputp=zeros(nr,nt);
        for c=1:nr
            routputx(c,:)=vx(xr(c),zr(c),:);
            routputz(c,:)=vz(xr(c),zr(c),:);
            for n=1:nt-1
                routputp(c,n)=(p(xr(c),zr(c),n)+p(xr(c),zr(c),n+1))/2;
            end
        end
end


L2Normzp=zeros(1,nr);L2Normxp=zeros(1,nr);
L2Normzs=zeros(1,nr);L2Normxs=zeros(1,nr);
L2Normpressure=zeros(1,nr);

% !! The norm of seismograms doesnt work in heterogeneous media!(the time window position is not correct)
if model==2
    disp('**!!!The norm of seismograms doesnt work in heterogeneous media!!!**(the time window position is not correct in the original codes)');
end




for c=1:nr
    norm1=max(abs([routputx(c,:),routputz(c,:)]));
    norm2=max(abs([as_ux(c,:),as_uz(c,:)]));
    routputx(c,:)=routputx(c,:)/norm1;
    routputz(c,:)=routputz(c,:)/norm1;
    as_ux(c,:)=as_ux(c,:)/norm2;
    as_uz(c,:)=as_uz(c,:)/norm2;



    norm1=max(abs([routputp(c,:)]));
    norm2=max(abs([as_p(c,:)]));
    routputp(c,:)=routputp(c,:)/norm1;
    as_p(c,:)=as_p(c,:)/norm2;


    L2Normxs(1,c)=norm(routputx(c,500:900)-as_ux(c,500:900),2);
    L2Normzs(1,c)=norm(routputz(c,500:900)-as_uz(c,500:900),2);
    L2Normxp(1,c)=norm(routputx(c,200:500)-as_ux(c,200:500),2);
    L2Normzp(1,c)=norm(routputz(c,200:500)-as_uz(c,200:500),2);
    % L2Normxp(1,c) = sqrt(norm(routputx(c,500:1000)-as_ux(c,500:1000),2))/sqrt(sum(as_ux(c,500:1000).^2));
    % L2Normzp(1,c) = sqrt(norm(routputz(c,500:1000)-as_uz(c,500:1000),2))/sqrt(sum(as_uz(c,500:1000).^2));
    % L2Normxs(1,c) = sqrt(norm(routputx(c,1:500)-as_ux(c,1:500),2))/sqrt(sum(as_ux(c,1:500).^2));
    % L2Normzs(1,c) = sqrt(norm(routputz(c,1:500)-as_uz(c,1:500),2))/sqrt(sum(as_uz(c,1:500).^2));
    L2Normpressure(1,c)=norm(routputp(c,:)-as_p(c,:),2);
end
as_outputx=as_ux(:,:);
as_outputz=as_uz(:,:);

kph=2*pi*fc*dh/vp0;
ksh=2*pi*fc*dh/vs0;
lengthkph=length(kph);
lengthksh=length(ksh);
cfdvp_c=zeros(lengthkph,nr);
cfdvs_c=zeros(lengthksh,nr);
PHI=0:0.01:2*pi-0.01; %phi is the angle between propagation direction and zaxis
if M==2
% L2norm_cfdvp_c=zeros(3,length(PHI));
% L2norm_cfdvs_c=zeros(3,length(PHI));
% L2_Normzp=zeros(3,nr);L2_Normxp=zeros(3,nr);
% L2_Normzs=zeros(3,nr);L2_Normxs=zeros(3,nr);
L2_Normpressure=zeros(3,nr);
end
% L2_Normzp(M/2,:)=L2Normzp(1,:);
% L2_Normxp(M/2,:)=L2Normxp(1,:);
% L2_Normzs(M/2,:)=L2Normzs(1,:);
% L2_Normxs(M/2,:)=L2Normxs(1,:);
L2_Normpressure(M/2,:)=L2Normpressure(1,:);
% for phi=[2*pi/nr:2*pi/nr:pi-2*pi/nr,pi+2*pi/nr:2*pi/nr:2*pi-2*pi/nr]  %phi is the angle between propagation direction and zaxis

for i=1:length(PHI)
    phi=PHI(i);
    if FD == 'CG'
        [cfdvp_c(:,i),~]=Z_Function_CG_homogeneous_cfd_c(phi,vp,dt,dh,M,kph,sourcetype,am,bm,max(max(gamma)));
        [~,cfdvs_c(:,i)]=Z_Function_CG_homogeneous_cfd_c(phi,vp,dt,dh,M,ksh,sourcetype,am,bm,max(max(gamma)));

    else
        [cfdvp_c(:,i),~]=Z_Function_SG_homogeneous_cfd_c(phi,vp,vs,dt,dh,M,kph,b_m);
        [~,cfdvs_c(:,i)]=Z_Function_SG_homogeneous_cfd_c(phi,vp,vs,dt,dh,M,ksh,b_m);
    end
    L2norm_cfdvp_c(M/2,i)=norm(1-cfdvp_c(:,i),2);
    L2norm_cfdvs_c(M/2,i)=norm(1-cfdvs_c(:,i),2);
end

if sourcetype==1
    if FD=='CG'
        char='acoustic ';
    else 
        char='';
    end
    %%  Seismogram acoustic
    figure(6);set(gcf, 'Position', [600, 600, 600, 400]);
    nexttile;
    for c=1:nr
        plot(t,routputp(c,:)+c,Color='k',linewidth=2.0,DisplayName='FD');hold on;
        xlim([0 nt*dt]);ylim([-2 nr+1]);grid on;hold on;
        if model==1
        plot(t,as_p(c,:)+c,'-',Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
        end
        title(['p seismograms,2M=',num2str(M)]);
        xlabel('time/s');
        ylabel('trace ID');
        if model==1
        lgd = legend('FD','GF','Location', 'southeast');
        else
            lgd = legend('FD','Location', 'southeast');
        end
        lgd.ItemTokenSize = [15, 18];
    end
    sg_str = sprintf('Normalized seismograms: %s %s-FD vs. analytical solutions' , ...
                  char,FD);
    if model==2
    sg_str = sprintf('Normalized seismograms of %s %s-FD' , ...
                  char,FD);
    end
    sgtitle(sg_str, 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Seismogram'],['p_seismogram_between_',FD,'FD_and_GF_2M=2_4_6.png']));
   
    
    %% Seismogram Norm acostic 
    figure(7);
    set(gcf, 'Position', [600, 600, 600, 400]);
    p1=plot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineStyle="none",LineWidth=2);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('L^2-Norm');
    % p1=polarplot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineWidth=1);
    hold on;
    legend('2M=2','2M=4','2M=6');
    if M==2
        p1.Marker='o';
    elseif M==4
        p1.Marker='+';
    elseif M==6
        p1.Marker='x';
    end
    title_str=sprintf('Seismogram error of %s%s-FD P-wave',char,FD);
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % ax = gca;
    % ax.ThetaDir = 'clockwise';
    % ax.ThetaZeroLocation = 'top';
    % saveas(gcf,'L2Norm_polarplot_between_CGFD_P_wave_pressure_seismograms_and_analytical_solutions.png');
    if M == 6
    ax_inset = axes('Position', [0.25, 0.2, 0.5, 0.2]);
    hold(ax_inset, 'on'); 
    grid(ax_inset, 'on');
    box(ax_inset, 'on');
    set(ax_inset, 'Color', 'w', 'FontSize', 10, 'FontName', 'Times New Roman');
    p1=plot(ax_inset, 0:2*pi/nr:2*pi, [L2_Normpressure(1,:),L2_Normpressure(1,1)], 'LineStyle', 'none', 'LineWidth', 2, 'DisplayName', ['2M=',num2str(M)]);
     p1.Marker='o';
    p1=plot(ax_inset, 0:2*pi/nr:2*pi, [L2_Normpressure(2,:),L2_Normpressure(2,1)], 'LineStyle', 'none', 'LineWidth', 2, 'DisplayName', ['2M=',num2str(M)]);
   p1.Marker='+';
    p1=plot(ax_inset, 0:2*pi/nr:2*pi, [L2_Normpressure(3,:),L2_Normpressure(3,1)], 'LineStyle', 'none', 'LineWidth', 2, 'DisplayName', ['2M=',num2str(M)]);
   p1.Marker='x';
    xlim(ax_inset, [0, 3.14]);
    maxlnormp = max(max(L2_Normpressure(2,:)),max(L2_Normpressure(3,:)));
    minlnormp = min(min(L2_Normpressure(2,:)),min(L2_Normpressure(3,:)));
    ylim(ax_inset, [minlnormp/1.05, maxlnormp*1.05]);title(ax_inset, 'Zoomed-in Detail', 'FontSize', 10);
    end


    
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Accuracy_analysis'],['L2Norm_plot_between_',FD,'FD_P_wave_pressure_seismograms_and_analytical_solutions.png']));



elseif sourcetype==2
    if FD=='CG'
        char='elastic ';
    else
        char='';
    end
    %% Elastic Seismogram x 
    figure(6);set(gcf, 'Position', [600, 600, 600, 400]);  
    % tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
    nexttile;hold off;
    % subplot(1,2,1);hold off;
    for c=1:nr
        plot(t,routputx(c,:)+c,Color='k',linewidth=2.0,DisplayName='FD');
        xlim([0 nt*dt]);ylim([-2.2 nr+1]);grid on;hold on;
        if model==1
        plot(t,as_outputx(c,:)+c,Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
        end
        if model==1
        lgd = legend('FD','GF','Location', 'southeast');
        else
            lgd = legend('FD','Location', 'southeast');
        end
        lgd.ItemTokenSize = [15, 18];
    end
    % title(['Normalized waveform ',char1,'_1,2M=',num2str(M)])
    title([char1,'_1 seismograms, 2M=',num2str(M)])
    xlabel('time/s');
    ylabel('trace ID')
    sg_str = sprintf('Normalized seismograms: %s%s-FD vs. analytical solutions' , ...
                  char,FD);
    if model==2
    sg_str = sprintf('Normalized seismograms of %s %s-FD' , ...
                  char,FD);
    end
    sgtitle(sg_str, 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Seismogram'],[FD,'x_seismograms_and_analytical_solutions','2M_246.png']));


 %% Elastic Seismogram z
    figure(7);set(gcf, 'Position', [600, 600, 600, 400]);  
    nexttile;hold off;
    % subplot(1,2,2);hold off;
    for c=1:nr
        plot(t,routputz(c,:)+c,Color='k',linewidth=2.0,DisplayName='FD');
        xlim([0 nt*dt]);ylim([-2.2 nr+1]);grid on;hold on;
        if model==1
        plot(t,as_outputz(c,:)+c,Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
        end
    end
    % title(['Normalized waveform ',char1,'_2,2M=',num2str(M)]);
    title([char1,'_2 seismograms, 2M=',num2str(M)])

    xlabel('time/s');
    ylabel('trace ID');
   if model==1
        lgd = legend('FD','GF','Location', 'southeast');
        else
            lgd = legend('FD','Location', 'southeast');
        end
    lgd.ItemTokenSize = [15, 18];

    sg_str = sprintf('Normalized seismograms: %s%s-FD vs. analytical solutions' , ...
                  char,FD);
    if model==2
    sg_str = sprintf('Normalized seismograms of %s %s-FD' , ...
                  char,FD);
    end
    sgtitle(sg_str, 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Seismogram'],[FD,'z_seismograms_and_analytical_solutions','2M_246.png']));


%%  Seismogram Norm elastic x
    figure(9);set(gcf, 'Position', [600, 600,600, 400]);
    subplot(1,2,1)
    p1=plot(0:2*pi/nr:2*pi,[L2Normxp,L2Normxp(1,1)],LineStyle="none",LineWidth=2);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('L^2-Norm');hold on;
    % p1=polarplot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineWidth=1);
    hold on;
    legend('2M=2','2M=4','2M=6');
    if M==2
        p1.Marker='o';
    elseif M==4
        p1.Marker='+';
    elseif M==6
        p1.Marker='x';
    end
     title_str=sprintf('P-wave error');
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % saveas(gcf,'L2Norm_plot_between_CGFD_P-wave_u_x_seismograms_and_analytical_solutions.png');
    % figure(7)
    subplot(1,2,2)
    % set(gcf, 'Position', [600, 600, 600, 600]);
    p2=plot(0:2*pi/nr:2*pi,[L2Normxs,L2Normxs(1,1)],LineStyle="none",LineWidth=2);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('L^2-Norm');hold on;
    % p1=polarplot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineWidth=1);
    hold on;
    legend('2M=2','2M=4','2M=6');
    if M==2
        p2.Marker='o';
    elseif M==4
        p2.Marker='+';
    elseif M==6
        p2.Marker='x';
    end
        
    title_str=sprintf('S-wave error');
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    sg_str=sprintf('%s_1 Seismogram error of %s%s-FD',char1,char,FD);
    sgtitle(sg_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % saveas(gcf,'L2Norm_plot_between_CGFD_S-wave_u_x_seismograms_and_analytical_solutions.png');
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Accuracy_analysis'],['L2Norm_plot_between_',FD,'FD_PandS-wave_',char1,'_1_seismograms_and_analytical_solutions.png']));
%%  Seismogram Norm elastic z

     figure(10);set(gcf, 'Position', [600, 600, 600,400]);
    subplot(1,2,1)
    p1=plot(0:2*pi/nr:2*pi,[L2Normzp,L2Normzp(1,1)],LineStyle="none",LineWidth=2);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('L^2-Norm');hold on;
    % p1=polarplot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineWidth=1);
    hold on;
    legend('2M=2','2M=4','2M=6');
    if M==2
        p1.Marker='o';
    elseif M==4
        p1.Marker='+';
    elseif M==6
        p1.Marker='x';
    end
    title_str=sprintf('P-wave error');
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % saveas(gcf,'L2Norm_plot_between_CGFD_P-wave_u_x_seismograms_and_analytical_solutions.png');
    % figure(7)
    subplot(1,2,2)
    % set(gcf, 'Position', [600, 600, 600, 600]);
    p2=plot(0:2*pi/nr:2*pi,[L2Normzs,L2Normzs(1,1)],LineStyle="none",LineWidth=2);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('L^2-Norm');hold on;
    % p1=polarplot(0:2*pi/nr:2*pi,[L2Normpressure,L2Normpressure(1,1)],LineWidth=1);
    hold on;
    legend('M=2','M=4','M=6');
    if M==2
        p2.Marker='o';
    elseif M==4
        p2.Marker='+';
    elseif M==6
        p2.Marker='x';
    end
    
    title_str=sprintf('S-wave error');
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    sg_str=sprintf('%s_2 Seismogram error of %s%s-FD',char1,char,FD);
    sgtitle(sg_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % saveas(gcf,'L2Norm_plot_between_CGFD_S-wave_u_x_seismograms_and_analytical_solutions.png');
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Accuracy_analysis'],['L2Norm_plot_between_',FD,'FD_PandS-wave_',char1,'_2_seismograms_and_analytical_solutions.png']));

end





%% Norm Dispersion curve P
figure(11);
set(gcf, 'Position', [0, 600, 600, 400]);
p3=plot(PHI,L2norm_cfdvp_c(M/2,:),LineStyle="-",LineWidth=1,DisplayName=['M=',num2str(M)]);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('$\|c_{fd}/c-1\|_2$',Interpreter='latex');
% p3=polarplot(0:0.01:2*pi,L2norm_cgcfdvp_c,LineWidth=1.0);hold on;
% legend('M=2','M=4','M=6');ylim([0 2.9e-3]);
legend show;
hold on;
% p3.LineStyle='none';
% if M==2
%     p3.Marker='o';
% elseif M==4
%     p3.Marker='+';
% elseif M==6
%     p3.Marker='x';
% end
% 

title_str=sprintf('Numerical dispersion error of %s%s-FD P-wave',char,FD);
title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
if M == 6
    ax_inset = axes('Position', [0.25, 0.2, 0.5, 0.2]);
    hold(ax_inset, 'on'); 
    grid(ax_inset, 'on');
    box(ax_inset, 'on');
    set(ax_inset, 'Color', 'w', 'FontSize', 10, 'FontName', 'Times New Roman');
    plot(ax_inset, PHI, L2norm_cfdvp_c(1,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['2M=',num2str(M)]);
    plot(ax_inset, PHI, L2norm_cfdvp_c(2,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['2M=',num2str(M)]);
    plot(ax_inset, PHI, L2norm_cfdvp_c(3,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['2M=',num2str(M)]);
    xlim(ax_inset, [0, 3.14]);
    maxlnormcfdvpc = max(max(L2norm_cfdvp_c(2,:)),max(L2norm_cfdvp_c(3,:)));
    minlnormcfdvpc = min(min(L2norm_cfdvp_c(2,:)),min(L2norm_cfdvp_c(3,:)));
    ylim(ax_inset, [minlnormcfdvpc/2, maxlnormcfdvpc*1.2]);title(ax_inset, 'Zoomed-in Detail', 'FontSize', 10);
end

% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% saveas(gcf,'L2Norm_polarplot_between_CGFD_P_wave_pressure_dispersion_relation_cfd_c_and_one.png');
saveas(gcf,fullfile([ResultFig,'/',FD,'_Accuracy_analysis'],['L2Norm_plot_between_',FD,'FD_',char,'P_wave_dispersion_relation_cfd_c_and_one.png']));

if sourcetype~= 1
    %% Norm dispersion curve elastic S
    figure(12);
    set(gcf, 'Position', [0, 600, 600, 400]);
    p3=plot(PHI,L2norm_cfdvs_c(M/2,:),LineStyle="-",LineWidth=1,DisplayName=['M=',num2str(M)]);grid on;axis tight;xlim([0 2*pi]);xlabel('θ/rad');ylabel('$\|c_{fd}/c-1\|_2$','Interpreter','latex');
    % p3=polarplot(0:0.01:2*pi,L2norm_cgcfdvp_c,LineWidth=1.0);hold on;
    % legend('M=2','M=4','M=6');ylim([0 2.9e-3]);
    legend show;
    hold on;
    % p3.LineStyle='none';
    % if M==2
    %     p3.Marker='o';
    % elseif M==4
    %     p3.Marker='+';
    % elseif M==6
    %     p3.Marker='x';
    % end
    title_str=sprintf('Numerical dispersion error of %s%s-FD S-wave',char,FD);
    title(title_str,'Interpreter','tex','FontSize', 12, 'FontWeight', 'bold');
    % ax = gca;
    % ax.ThetaDir = 'clockwise';
    % ax.ThetaZeroLocation = 'top';
    % saveas(gcf,'L2Norm_polarplot_between_CGFD_P_wave_pressure_dispersion_relation_cfd_c_and_one.png');
    if M == 6
        if strcmp(FD, 'CG') && strcmp(char, 'elastic ')
            ax_inset = axes('Position', [0.5, 0.4, 0.4, 0.2]);
        else
            ax_inset = axes('Position', [0.25, 0.2, 0.5, 0.2]);
        end
    hold(ax_inset, 'on'); 
    grid(ax_inset, 'on');
    box(ax_inset, 'on');
    set(ax_inset, 'Color', 'w', 'FontSize', 10, 'FontName', 'Times New Roman');
    plot(ax_inset, PHI, L2norm_cfdvs_c(1,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['M=',num2str(M)]);
    plot(ax_inset, PHI, L2norm_cfdvs_c(2,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['M=',num2str(M)]);
    plot(ax_inset, PHI, L2norm_cfdvs_c(3,:), 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', ['M=',num2str(M)]);
    xlim(ax_inset, [0, 3.14]);
    maxlnormcfdvsc = max(max(L2norm_cfdvs_c(2,:)),max(L2norm_cfdvs_c(3,:)));
    minlnormcfdvsc = min(min(L2norm_cfdvs_c(2,:)),min(L2norm_cfdvs_c(3,:)));
    ylim(ax_inset, [minlnormcfdvsc/2, maxlnormcfdvsc*1.2]);title(ax_inset, 'Zoomed-in Detail', 'FontSize', 10);
    end
    saveas(gcf,fullfile([ResultFig,'/',FD,'_Accuracy_analysis'],['L2Norm_plot_between_',FD,'FD_',char,'S_wave_dispersion_relation_cfd_c_and_one.png']));
end





% subplot(2,2,2);
% p2=polarplot(0:2*pi/nr:2*pi,[L2Normxs,L2Normxs(1,1)],LineWidth=1);hold on;legend('M=2','M=4','M=6');
% p2.LineStyle='none';
% if M==2
%     p2.Marker='o';
% elseif M==4
%     p2.Marker='+';
% elseif M==6
%     p2.Marker='x';
% end
% title('L2Norm polarplot between CGFD u_x S-wave seismograms and analytical solutions');
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';


% subplot(2,2,4);
% p4=polarplot(0:2*pi/nr:2*pi,[L2Normzs,L2Normzs(1,1)],LineWidth=1);hold on;legend('M=2','M=4','M=6');
% p4.LineStyle='none';
% if M==2
%     p4.Marker='o';
% elseif M==4
%     p4.Marker='+';
% elseif M==6
%     p4.Marker='x';
% end
% title('L2Norm polarplot between CGFD u_y S-wave seismograms and analytical solutions');
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';


% figure(4);
% set(gcf, 'Position', [650, 0, 1300, 1200]);
% plot(0:2*pi/nr:2*pi-2*pi/nr,L2Normz(1,:),'*','LineWidth',2);grid on;hold on;legend('M=2','M=4','M=6');
% title('L2Norm plot between CGFD u_z seismograms and analytical solutions');
% ylim([0 2.2]);
% xlabel('Clockwise angle from top/rad');
% ylabel('L2Norm');
% saveas(gca,'CGL2Normpolt.png');

%
%  out=ux(:,:,nt);
%  figure(2);
%  set(gcf, 'Position', [0, -100, 1080, 980]);
%  subplot(2,2,1);
%  imagesc(out');
% % set(gcf, 'Position', [0, 0, 1300, 1200]);
%  hold on;
%  for c=1:nr
%      plot(xr(c),zr(c),'*',Color='b');
%  end
% plot(xs,zs,'*',Color='r') ;
% axis tight;
% colormap(jet);
% colorbar;
% out=uz(:,:,nt);
%  figure(2);
%  set(gcf, 'Position', [0, -100, 1080, 980]);
%  subplot(2,2,3);
%  imagesc(out');
% % set(gcf, 'Position', [0, 0, 1300, 1200]);
%  hold on;
%  for c=1:nr
%      plot(xr(c),zr(c),'*',Color='b');
%  end
% plot(xs,zs,'*',Color='r') ;
% axis tight;
% colormap(jet);
% colorbar;
% figure(2);
% subplot(2,2,3);
% for c=1:nrp
%     plot(t,routputx(c,:)-as_outputx(c,:));hold on;
%
% end
% for c=1:nr
%     plot(t,routputz(c,:)-as_outputz(c,:)+1);hold on;
%
% end
% xlim([0 nt*dt]);ylim([-0.5 1.5]);
% grid on;
% hold on;
% L1norm=norm(routputz(c,:)-as_outputz(c,:),1);
% L2norm=norm(routputz(c,:)-as_outputz(c,:),2);
% disp(['M=',num2str(M)]);
% disp(['L1norm= ',num2str(L1norm)]);
% disp(['L2norm= ',num2str(L2norm)]);









