function [vx,vy,t]=GF2Del(fc,vp0,vs0,rho,dt,T,x,y,FD)
%fc center freq 
%vp0 velocity p wave
%vs0 velocity s wave\
%rho density
%dt time interval
%T total time 
%x y
n2ft=T/dt;%n2ft number of time samples(no 0 time here)

t=(1:n2ft)*dt; %t time series(no 0 here)
fnyq=1/(2*dt); %fnyq nyquist freq

%here define the Ricker wave
%%tsour=1/fc;%td=1/fc td domiant period



%%tau=pi*(t(1:(2*round(tsour/dt)))-tsour)/tsour;%this one could not work
%because it is renewed later


%tau=pi(t-td)fc
%t-delay means time domain moves left %renew tau here 
delay=1.0/fc; shift=0.0/fc; ts=t-delay-shift; tau=pi*fc*ts;
s=(1.0-2.0*tau.*tau).*exp(-tau.*tau);%(1-2tau^2)e^{-2tau^2}
if FD == 'CG'
    s=tau.*exp(-tau.*tau); 
end
%%t0=delay;%this one is not used

S=fftshift(fft(s,n2ft));%change to freq domain

f=fnyq*(-n2ft/2:n2ft/2-1)/(n2ft/2); %freq series
w=2*pi*f;%change to angle freq
w(n2ft/2+1)=w(n2ft/2);

vp=zeros(n2ft)+vp0;%vp in time series

vs=zeros(n2ft)+vs0;%vs in time series

F=1.0;

% correct for staggered FD position of vx shifted by +dh/2 in x
% and -dh/2 in z
% ony for x-component u1
%dh=1.0;

% x2=x;
% y2=y-dh/2;
% x1=x-dh/2;
% y1=y-dh/2;
% x1=x-dh/2;
% y1=y+dh/2;
x1=x;x2=x;
y1=y;y2=y;


r1=sqrt(x1^2+y1^2);     % Distance to the receiver for u1
r2=sqrt(x2^2+y2^2);     % Distance to the receiver for u2


%% Green's functions (CORRECTED, as appear in Carcione (2002))
for j=1:length(w)
    if w(j)>0
        G1=-1i*pi/2*(1/vp(j)^2*besselh(0,2,(w(j)*r1/(vp(j))))+ ...
            1./(w(j)*r1*vs(j))*besselh(1,2,(w(j)*r1/(vs(j))))- ...
            1./(w(j)*r1*vp(j))*besselh(1,2,(w(j)*r1/(vp(j)))));
        G2=1i*pi/2*(1/vs(j)^2*besselh(0,2,(w(j)*r2/(vs(j))))- ...
            1./(w(j)*r2*vs(j))*besselh(1,2,(w(j)*r2/(vs(j))))+ ...
            1./(w(j)*r2*vp(j))*besselh(1,2,(w(j)*r2/(vp(j)))));
        u1(j)=F/(2*pi*rho)*(x1*y1/r1^2)*(G1+G2);
        u2(j)=F/(2*pi*rho)*(1/r2^2)*(y2^2*G1-x2^2*G2);
    end
end



%Ricker type source
%S=sqrt(pi)*w.^2/(4*(pi*fc)^3).*exp(-i*w*t0).*exp(-w.^2/(4*(pi*fc)^2));


%% Construction of the w-space displacements
PHI_DX=u1.*S;
PHI_DY=u2.*S;

%% Construction of the w-space velocities
PHI_VX=PHI_DX.*i.*w;
PHI_VY=PHI_DY.*i.*w;

%% FFT of the PHI functions to obtain solution
Sol_x=real((ifft(fftshift(PHI_DX))));  
Sol_y=real((ifft(fftshift(PHI_DY))));

%% Same for the velocities 
Vel_x=real((ifft(fftshift(PHI_VX))))/dt;  
Vel_y=real((ifft(fftshift(PHI_VY))))/dt;


% vx=Sol_x;
% vy=Sol_y;
vx=Vel_x;
vy=Vel_y;

if FD=='SG'
    vx=Sol_x;
    vy=Sol_y;
end




%%%%%%%%%%%%%%%% 
%% REFERENCES %%
%%%%%%%%%%%%%%%%
% (1) "WAVE PROPAGATION SIMULATION IN A LINEAR VISCOELASTIC MEDIUM", Carcione et al. (1988)
%
% (2) "WAVE FIELDS IN REAL MEDIA: WAVE PROPAGSTION IN ANISOTROPIC, ANELASTIC
% AND POROUS MEDIA", Carcione (2002)
%
% (3) "THE FINITE DIFFERENCE METHOD FOR SEISMOLOGISTS: AN INTRODUCTION",
% Moczo et al. (2005)
%
% (4) "INCORPORATION OF ATTENUATION INTO TIME-DOMAIN COMPUTATIONS OF
% SEISMIC WAVE FIELDS", Emmerich & Korn (1987)

