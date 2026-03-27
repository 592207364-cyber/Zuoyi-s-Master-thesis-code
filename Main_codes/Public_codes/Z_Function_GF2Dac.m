function [p,t]=GF2Dac(fc,c,R,dt,T)
%fc center freq 
%c velocity. k wavenumber w =2pif angel freq 
%dt time intervial
%T total time





n2ft=T/dt;%n2ft number of time samples(no 0 time here)

t=(1:n2ft)*dt; %t time series(no 0 here)
fnyq=1/(2*dt); %fnyq nyquist freq

%here define the Ricker wave
%%tsour=1/fc;%td=1/fc td domiant period

%%tau=pi*(t(1:(2*round(tsour/dt)))-tsour)/tsour; %this one could not work
%because it is renewed later

%tau=pi(t-td)fc
%t-delay means time domain moves left %renew tau here 
shift=0.0/fc;pt0=Z_Function_Ricker_wave(t,fc,shift);
% delay=1.0/fc; shift=0.0/fc; ts=t-delay-shift; tau=pi*fc*ts;
% pt0=(1.0-2.0*tau.*tau).*exp(-tau.*tau);%(1-2tau^2)e^{-2tau^2}

PT0=fftshift(fft(pt0,n2ft));%change to freq domain

f=fnyq*(-n2ft/2:n2ft/2-1)/(n2ft/2); 
%give a freq series//but I cannot understand this now
% it could be to avoid 0 freq? 

w=2*pi*f;%change to angle freq series %w cannot be 0
w(n2ft/2+1)=w(n2ft/2);
k=w./c;%change to wavenumber 
 
I2D=sqrt(2.0*pi./(k*R));%G2d=sqrt(2pic/wr) e^-ikr e^-pi/4
I2D=I2D.*exp(-1i*k*R-pi/4.0);%solution of GF2d 

P2D=I2D.*PT0; %freq multiply=time domain convolution

p2D=real(ifft(ifftshift(P2D),n2ft));%ifft
p=p2D;
%so.,this function is to give a analytic solution of GF2D Ricker wave
