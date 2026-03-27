function [sgcfdvp_c,sgcfdvs_c]=Z_Function_SG_homogeneous_cfd_c(phi,vp,vs,dt,dh,M,kh,b_m)
% kh=0:0.01:pi;
lengthkh=length(kh);
rp=vp*dt/dh; rp=max(max(rp));
rs=vs*dt/dh; rs=max(max(rs));
sgcfdvp_c=zeros(1,lengthkh);
sgcfdvs_c=zeros(1,lengthkh);
for j=1:lengthkh
    sum1=0;sum2=0;
    for m=1:M/2
        sum1=sum1+b_m(m)*sin(kh(j)*cos(phi)*(m-0.5));
        sum2=sum2+b_m(m)*sin(kh(j)*sin(phi)*(m-0.5));
    end
    sgcfdvp_c(j)=2/kh(j)/rp*asin(rp*sqrt(sum1^2+sum2^2));
    sgcfdvs_c(j)=2/kh(j)/rs*asin(rs*sqrt(sum1^2+sum2^2));
end
