i=1;
  figure(2);
   set(gcf, 'Position', [0, 0, 700, sourcetype*350]);
 tiledlayout(sourcetype,2,"TileSpacing","tight",Padding="tight");
for n=1:nt
    for c=1:nr
        switch FD
            case 'CG'
                routputx(c,n)=ux(xr(c),zr(c),n)+c;
                routputx(c,n)=uz(xr(c),zr(c),n)+c;
                normx=max(max(abs([ux(:,:,n)])));
                normz=max(max(abs([uz(:,:,n)])));
            case 'SG'
                routputx(c,n)=vx(xr(c),zr(c),n)+c;
                routputz(c,n)=vz(xr(c),zr(c),n)+c;
                normx=max(max(abs([vx(:,:,n)])));
                normz=max(max(abs([vz(:,:,n)])));
        end
    end
    %% image
    if n==500 || n==666


        %% image snapshot p or v_x or u_x
        figure(2);
         nexttile;
        %subplot(2,2,i)
        i=i+1;
        if sourcetype==1
            normx=max(max(abs([p(:,:,n)])));
            outx=p(:,:,n)/normx;
            char1='Acoustic';
            char2='p';
        else
            char1='Elastic';
            switch FD
                case 'CG'
                    outx=ux(:,:,n)/normx;
                   
                    char2='u_1';
                case 'SG'
                    outx=vx(:,:,n)/normx; 
                    char2='v_1';

            end
        end
        imagesc(outx');hold on;
        if model==2
           % line([xs+squarex1,xs+squarex2,xs+squarex2,xs+squarex1,xs+squarex1],[zs+squarez1,zs+squarez1,zs+squarez2,zs+squarez2,zs+squarez1],'Color','r','LineWidth',2);
            line([xs+squarex1,xs+squarex2],[zs+squarez1,zs+squarez1],'Color','k','LineWidth',1);
        end
        % set(gcf, 'Position', [0, 0, 1000, 300]);%first 2 value is the location at the monitor.
        plot(xr, zr, 'b*', 'MarkerSize', 6); 
        plot(xs, zs, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5); 
        plot(xs,zs,'*',Color='r') ;hold on;
        %to check which receiver is the first one and the last one
        %plot(xr(1),zr(1),'*',Color='g');plot(xr(nr),zr(nr),'*',Color='k');
        title([char2,' snapshot at t = ',num2str(n*dt),' s']);
            % title(['t = ',num2str(n*dt),' s, 2M=',num2str(M)]);
        xlabel('$x_1$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$x_2$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
        pbaspect([1 1 1]);
        axis tight;
        cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
        colormap(gca, cmap);
        if i==3 %% let the last figure have colorbar
            colorbar;
        end
        drawnow;
        clim([-1 1]);
        if sourcetype==1
        saveas(gcf,fullfile([ResultFig,'/',FD,'_Snapshots'],[FD,'_',char1,'_snapshots_with_2M=',num2str(M),'.png']));
        end
    end

    %% This part is for the PPT FIG
        % if mod(n,666)==0
        %     figure(9);
        %     nexttile;
        %     imagesc(outx');hold on;
        % 
        %     set(gcf, 'Position', [0, 0, 1000, 300]);%first 2 value is the location at the monitor.
        %     plot(xr, zr, 'b*', 'MarkerSize', 6); 
        %     plot(xs, zs, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5); 
        %     plot(xs,zs,'*',Color='r') ;hold on;
        %     %to check which receiver is the first one and the last one
        %     %plot(xr(1),zr(1),'*',Color='g');plot(xr(nr),zr(nr),'*',Color='k');
        %     title([FD,char1,' snapshot at t = ',num2str(n*dt),' s',', 2M = ',num2str(M)]);
        %     xlabel('$x_1$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
        %     ylabel('$x_2$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
        %     pbaspect([1 1 1]);
        %      axis tight;
        %     cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
        %     ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
        %      colormap(gca, cmap);
        %     if M==6 %% let the last figure have colorbar
        %         colorbar;
        %     end
        %     drawnow;
        %     clim([-1 1])
        %     saveas(gcf,fullfile([ResultFig,'/',FD,'_Snapshots'],[FD,'_',char1,'_snapshots_with_2M=246','.png']));
        % end
    % This part is for the PPT FIG(end)
end

if sourcetype~=1
    switch FD
        case 'CG'
            char2='u_2';
            % char3='u_1';
        case 'SG'
            char2='v_2';
            % char3='v_1';
    end
    i=1;
    % figure(3);
    % tiledlayout(1,3,"TileSpacing","tight",Padding="tight")
    
    for n=1:nt
        for c=1:nr
            switch FD
                case 'CG'
                    routputz(c,n)=uz(xr(c),zr(c),n)+c;
                    normz=max(max(abs([uz(:,:,n)])));
                case 'SG'
                    routputz(c,n)=vz(xr(c),zr(c),n)*20+c;
                    normz=max(max(abs([vz(:,:,n)])));
            end
        end
        if n==500 || n==666
            % figure(3);
            % subplot(1,3,i);
            % nexttile;
            figure(2)
            % subplot(2,2,i+2)
            nexttile;
             i=i+1;
            switch FD
                case 'CG'
                    outz=uz(:,:,n)/normz;
                case 'SG'
                    outz=vz(:,:,n)/normz;
            end
            imagesc(outz');hold on;
            % set(gcf, 'Position', [0, 0, 1000, 300]);%first 2 value is the location at the monitor.
            plot(xr, zr, 'b*', 'MarkerSize', 6); 
            plot(xs, zs, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5); 
            plot(xs,zs,'*',Color='r') ;
            if model==2
           % line([xs+squarex1,xs+squarex2,xs+squarex2,xs+squarex1,xs+squarex1],[zs+squarez1,zs+squarez1,zs+squarez2,zs+squarez2,zs+squarez1],'Color','r','LineWidth',2);
            line([xs+squarex1,xs+squarex2],[zs+squarez1,zs+squarez1],'Color','k','LineWidth',1);
            end
            title([char2,' snapshot at t = ',num2str(n*dt),' s']);
            % title(['t = ',num2str(n*dt),' s, 2M=',num2str(M)]);
            xlabel('$x_1$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('$x_2$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
            pbaspect([1 1 1]);
            axis tight;
            cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
            colormap(gca, cmap);
            if i==3
                colorbar;
            end
            drawnow;
            hold on;
            clim([-1 1]);
            sg_str = sprintf('%s %s-FD snapshots, 2M=%d' , ...
                  char1,FD,M);
            % sg_str = sprintf('%s %s-FD %s snapshots' , ...
                  % char1,FD,char3);
             sgtitle(sg_str, 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
    
            saveas(gcf,fullfile([ResultFig,'/',FD,'_Snapshots'],[FD,'_',char2,'_snapshots_with_2M=',num2str(M),'.png']));
        end
    end
end
