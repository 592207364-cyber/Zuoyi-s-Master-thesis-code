function s=ricker(t,fc,shift)

    ns=length(t);
    s=zeros(ns);
    delay=1.0/fc; ts=t-delay-shift; tau=pi*fc*ts;


     s= (1-2.*tau.^2).*exp(-tau.^2);
     %s=tau.*exp(-tau.*tau);
end
