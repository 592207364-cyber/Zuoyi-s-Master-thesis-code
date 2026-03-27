disp('Notice: The stability limits are obtained from the code - Grid Scan');
load("rmaxlimit.mat")
gammamax=max(max(gamma));
gammaindex=gammamax*10;%gamma M kh
if max(max(F)) <= rmaxlimit(gammaindex,M/2)
    disp('r<=rmax,stable');
else
    disp('r>rmax,UNSTABLE');
end