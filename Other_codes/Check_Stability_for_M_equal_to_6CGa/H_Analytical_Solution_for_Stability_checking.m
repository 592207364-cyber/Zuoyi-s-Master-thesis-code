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

as_outputx=as_ux(:,:);
as_outputz=as_uz(:,:);

kh=2*pi*fc*dh/vp0;
lengthkh=length(kh);
cfdvp_c=zeros(lengthkh,nr);
cfdvs_c=zeros(lengthkh,nr);
PHI=0:0.01:2*pi-0.01; %phi is the angle between propagation direction and zaxis
L2norm_cfdvp_c=zeros(1,length(PHI));
L2norm_cfdvs_c=zeros(1,length(PHI));
% for phi=[2*pi/nr:2*pi/nr:pi-2*pi/nr,pi+2*pi/nr:2*pi/nr:2*pi-2*pi/nr]  %phi is the angle between propagation direction and zaxis



if sourcetype==1
    char='acoustic';
    figure(3);
    set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 500]);
    tlayout = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
    nexttile;
    plot(t, routputp' + (1:nr), 'k', 'LineWidth', 1.5);
    grid on; hold on;
    xlim([0, nt*dt]); ylim([0, nr+1]);
    title(['Seismogram of p(unnormalized), 2M=', num2str(M), ', v_p=',num2str(vp0),' m/s'], 'FontSize', 11);
    xlabel('Time (s)'); ylabel('Trace ID');
    lgd = legend('FD', 'Location', 'southwest');
    lgd.ItemTokenSize = [15, 10];
    nexttile;
    plot(t, routputp'*10 + (1:nr), 'k', 'LineWidth', 1.5);
    grid on; hold on;
    xlim([nt*dt*0.75, nt*dt]); ylim([0, nr+1]);
    title(['Zoomed-in (Ampl. x10) p, 2M=', num2str(M),', v_p=',num2str(vp0),' m/s'], 'FontSize', 11);
    xlabel('Time (s)'); ylabel('Trace ID');
    lgd = legend('FD', 'Location', 'southwest');
    lgd.ItemTokenSize = [15, 10];
    exportgraphics(gcf, fullfile([FD,'_p_seismograms_verification.png']), 'Resolution', 300);

elseif sourcetype==2
    char='elastic';
    %Seismogram
    figure(M/2+5);
    tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
    nexttile;hold off;
    % subplot(1,2,1);hold off;
    for c=1:nr
        plot(t,routputx(c,:)+c,Color='k',linewidth=2.0,DisplayName='FD');
        xlim([0 nt*dt]);ylim([-1 nr+1]);grid on;hold on;
        % plot(t,as_outputx(c,:)+c,Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
    end
    title(['Normalized waveform ',char1,'_x,2M=',num2str(M)])
    xlabel('time/s');
    ylabel('trace ID')
    l=legend('FD','GF','Location','southwest');
    l.ItemTokenSize=[20,10];
    nexttile;hold off;
    % subplot(1,2,2);hold off;
    for c=1:nr
        plot(t,routputz(c,:)+c,Color='k',linewidth=2.0,DisplayName='FD');
        xlim([0 nt*dt]);ylim([-1 nr+1]);grid on;hold on;
        plot(t,as_outputz(c,:)+c,Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
    end
    title(['Normalized waveform ',char1,'_y,2M=',num2str(M)]);
    xlabel('time/s');
    ylabel('trace ID');
    l=legend('FD','GF','Location','southwest');
    l.ItemTokenSize=[20,10];
    % saveas(gcf,fullfile([ResultFig,'/',FD,'_Seismogram'],[FD,'_seismograms_and_analytical_solutions','2M_',num2str(M),'.png']));
    figure(12);
    nexttile;hold off;
    % subplot(1,2,1);hold off;
    for c=1:nr
        plot(t,routputx(c,:)*5+c,Color='k',linewidth=2.0,DisplayName='FD');
        xlim([0 nt*dt]);ylim([-1 nr+1]);grid on;hold on;
        % plot(t,as_outputx(c,:)+c,Color='r',linewidth=2.0,DisplayName='GF');%*1e12*0.4
    end
    title(['Unnormalized waveform ',char1,'_x,2M=',num2str(M)])
    xlabel('time/s');
    ylabel('trace ID')
    l=legend('FD','Location','southwest');
    l.ItemTokenSize=[20,10];
    saveas(gcf,fullfile([FD,'_',char1,'_x_snapshots_and_seismograms','.png']));
                


    if sourcetype==1
            norm=max(max(abs([p(:,:,n)])));
            out=p(:,:,n)/norm;
            char1='pressure';
        else
            switch FD
                case 'CG'
                    out=ux(:,:,n)/norm;
                    char1='u_x';
                case 'SG'
                    out=vx(:,:,n)/norm;
                    char1='v_x';
            end
        end
   
    


    





end