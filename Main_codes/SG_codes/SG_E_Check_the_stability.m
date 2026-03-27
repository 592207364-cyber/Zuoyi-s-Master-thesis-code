sumabsb_m=0;
for m=1:M/2
    sumabsb_m=sumabsb_m+b_m(m);
end
rmax=1/(sqrt(2)*sumabsb_m);
if F <= rmax
    disp([FD,',Stable'])
else
    disp([FD,',*UNSTABLE*'])
end
