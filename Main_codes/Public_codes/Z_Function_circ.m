function [xr,zr,phi]=Z_Function_circ(nr,r)
dphi=2*pi/nr;%angle split by receivers
n=1;%step
xr=zeros(1,nr);
zr=zeros(1,nr);
phi=zeros(1,nr);
for p=0:dphi:2*pi-dphi
    phi(n)=p;%angle
    xr(n)=r*sin(p);%x coordinate    
    zr(n)=r*cos(p);%z coordinate
    n=n+1;%next step
end
end