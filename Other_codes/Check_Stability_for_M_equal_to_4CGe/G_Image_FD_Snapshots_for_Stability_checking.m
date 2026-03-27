i=1;
figure(2);
tiledlayout(1,3,"TileSpacing","tight",Padding="tight")

for n=1:nt
    for c=1:nr
        switch FD
            case 'CG'
                routputx(c,n)=ux(xr(c),zr(c),n)+c;
                norm=1;
            case 'SG'
                routputx(c,n)=vx(xr(c),zr(c),n)*20+c;
                norm=1;
        end
    end
    if sourcetype==1
            norm=max(max(abs([p(:,:,n)])));
            out=p(:,:,n)/norm;
            char1='pressure';
        else
            switch FD
                case 'CG'
                    out=ux(:,:,n)/norm;
                    char1='u_1';
                case 'SG'
                    out=vx(:,:,n)/norm;
                    char1='v_1';
            end
        end
    %% image
    if mod(n,80)==0
        figure(2);
        nexttile;
        i=i+1;
        
        imagesc(out');hold on;
        set(gcf, 'Position', [0, 0, 1000, 300]);%first 2 value is the location at the monitor.
        for c=1:nr %plot receiver points
            plot(xr(c),zr(c),'*',Color='b');  hold on;
        end
        plot(xs,zs,'*',Color='r') ;hold on;
        %to check which receiver is the first one and the last one
        %plot(xr(1),zr(1),'*',Color='g');plot(xr(nr),zr(nr),'*',Color='k');
        title([FD,char1,' snapshot at t = ',num2str(n*dt),' s',', 2M = ',num2str(M)]);
        xlabel('X(m)',LineWidth=2);ylabel('Z(m)',LineWidth=2);
        axis tight;
        colormap(jet);
        if i==4 %% let the last figure have colorbar
            colorbar;
        end
        drawnow;
        clim([-1 1]);
        
        % saveas(gcf,fullfile([FD,'_Snapshots'],[FD,'_',char1,'_snapshots_with_2M=',num2str(M),'.png']));
    end

    %% This part is for the PPT FIG
     if mod(n,80)==0 && n~=80
    figure(12);
    set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 400]); 
    
    nexttile;
    imagesc(out'); hold on;

    plot(xr, zr, 'b*', 'MarkerSize', 6);  
    plot(xs, zs, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5); 

    title(sprintf('t = %.3f s', n*dt), 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');

    xlabel('$x_1$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$x_2$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    
    pbaspect([1 1 1]); 
    axis tight;

    cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
    colormap(gca, cmap);
    colorbar;
    clim([-0.2 0.2]); 

    sg_str = sprintf('Elastic %s-FD %s snapshot (unnormalized), 2M = %d, v_p=%.1f m/s, v_s=%.1f m/s', ...
                  FD, char1, M, vp0, vs0);
    sgtitle(sg_str, 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow;
    
    img_name = fullfile([FD, '_', char1, '_snapshots.png']);
    exportgraphics(gcf, img_name, 'Resolution', 600);
end
    % This part is for the PPT FIG(end)

          





end
