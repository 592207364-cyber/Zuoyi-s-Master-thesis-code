%% Accuracy order
% M=2;% accuracy order is M.
%% Define x-z grid parameters
nx=260;     % Model size in grid points
nz=260;     %#
dh=1.0;     %m
dz=dh;
dx=dh;
x=dx:dx:nx*dx;%x coordinate
z=dz:dz:nx*dz;%z coordinate
%% Define time parameters
dt=4.0e-4;%s
T=0.1; %s
nt=T/dt;%#
t=dt:dt:T; %(dt~nt*dt)   time coordinate
%% Source
xs=round(nx/2);			% source position in grid points (mid grid)
zs=round(nz/2);        %#
%% Receivers
nr=24; radius=50;
% calculate source and receiver coordinates
[xr,zr,phir]=Z_Function_circ(nr,radius);
% Discrete FD source positions
xr=round(xs+xr);
zr=round(zs+zr);%change xr zr from the raletive position of source to grid coodinate
%% receivers output
routputx=zeros(nr,nt);routputz=zeros(nr,nt);
%% Define model
rho0=2000.0;
vp=zeros(nx,nz)+vp0;             %m/s P wave velocity
vs=zeros(nx,nz)+vs0;               %m/s S wave velocity
rho=zeros(nx,nz)+rho0;             %density kg/m³
alpha=zeros(nx,nz); alpha=alpha+vp;  %P wave velocity
beta=zeros(nx,nz); beta=beta+vs;%S wave velocity field
%% Define Initial parameters for the equation
F=alpha.*dt./dh; %CFL Number% alpha beta given % dt dh given % dz dx given
mu=beta.^2.*rho;
if sourcetype==1
    mu=zeros(nx,nz);
end
lambda=alpha.^2.*rho-2.*mu;%Lamé parameter λ μ
fc=80;   %Hz
switch FD
    case 'CG'
        %% Define Ricker wave

        rickerf=Z_Function_Ricker_wave(t,fc,1.0*dt); % Force source
        rickerp=Z_Function_Ricker_wave(t,fc,2.0*dt); % pressure source
        % figure(1)
        % plot(t,rickerp);grid on;hold on;
        % plot(t,rickerf);legend('pressure source wavelet','force source wavelet');
        % title(['Ricker wavelet waveform,f_c=',num2str(fc),'Hz']);
        % xlabel('time/s')
        % ylabel('Amplitude')
        % saveas(gcf,'Ricker_wavelet_waveform.png');
        % hold off;
        %% Other parameters
        gamma=beta./alpha;%ratio vs/vp
        ux=zeros(nx,nz,nt); %u_x displacement_x
        uz=zeros(nx,nz,nt); %u_z displacement_z
        pxx=zeros(nx,nz,nt);
        pzz=zeros(nx,nz,nt);
        p=zeros(nx,nz,nt);
    case 'SG'
        %% Define Ricker wave
        rickerf=Z_Function_integrated_Ricker_wave(t,fc,1.0*dt); % Force source
        rickerp=Z_Function_integrated_Ricker_wave(t,fc,2.0*dt); % Pressure source
        % figure(1)
        % plot(t,rickerp);grid on;hold on;
        % plot(t,rickerf);legend('pressure source wave','force source wave');
        % title(['Ricker wavelet waveform,f_c=',num2str(fc),'Hz']);
        % xlabel('time/s')
        % ylabel('Amplitude')
        % saveas(gcf,'integrated_Ricker_wavelet_waveform.png');
        % hold off;
        lambda2mu=alpha.^2.*rho;
        lambdadt_dh=lambda.*dt./dh;%L4 Page11/34 first fomula lamdadt/h
        dt_rhoh=dt./dh./rho;%L4 Page 11/34 second fomula dt/rhoh
        % Initialize wavefield matrix
        %%2D staggered grid
        %%see lecture PDF page 21/34 shows all the parameter we need
        vx=zeros(nx,nz,nt);  %%V_x
        vz=zeros(nx,nz,nt);  %%V_z
        sxx=zeros(nx,nz,nt); %%sigma xx
        sxz=zeros(nx,nz,nt); %%sigma xz
        szz=zeros(nx,nz,nt); %%sigma zz
        p=zeros(nx,nz,nt);    %%p(i,j)=sxx(i,j)+szz(i,j);
end