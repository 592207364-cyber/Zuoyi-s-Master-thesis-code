function [st1,st2]=Z_Function_CG_ElasticFD_Stability_rmax(M,gamma0,dh,k,am,bm)
dtheta=0.01;
% theta=0:dtheta:2*pi;
theta=0:dtheta:pi;
lengththeta=length(theta);
sum1=zeros(lengththeta,1);sum2=zeros(lengththeta,1);sum3=zeros(lengththeta,1);
for i=1:lengththeta%1 ---theta=0 here
    for m=1:M/2
        sum1(i)=sum1(i)+am(m)*(sin(m*dh*k*sin((i-1)*dtheta)/2))^2;
        sum2(i)=sum2(i)+am(m)*(sin(m*dh*k*cos((i-1)*dtheta)/2))^2;
        index=0;
        for l=1:M/2
            index=index+bm(l)*sin(m*k*dh*sin((i-1)*dtheta))*sin(l*k*dh*cos((i-1)*dtheta));
        end
        sum3(i)=sum3(i)+bm(m)*index;
    end
end
% A=rp^2.*sum1+rp^2*gamma0^2.*sum2;
A0=sum1+gamma0^2.*sum2;
% D=rp^2.*sum2+rp^2*gamma0^2.*sum1;
D0=sum2+gamma0^2.*sum1;
% B=rp^2*(1-gamma0^2)/4*sum3;
B0=(1-gamma0^2)/4*sum3;
% S1=0.5*(A+D)+0.5*sqrt((A-D).^2+4*B.^2);
% S2=0.5*(A+D)-0.5*sqrt((A-D).^2+4*B.^2);
%st is rpmax formula
st1=A0+D0+sqrt((A0-D0).^2+4*B0.^2);st1=2./st1;st1=sqrt(st1);
st2=A0+D0-sqrt((A0-D0).^2+4*B0.^2);st2=2./st2;st2=sqrt(st2);
end