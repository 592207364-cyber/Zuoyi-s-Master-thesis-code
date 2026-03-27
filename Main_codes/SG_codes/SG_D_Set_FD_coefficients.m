%% FD coefficient bm=beta_m
A=zeros(M/2);
for l=1:M/2
    for k=1:M/2
        A(l,k)=(2*k-1)^(2*l-1);%(2m-1)^{2k-1}
    end
end
c=zeros(1,M/2);
c(1)=1;
c=c';
b_m=inv(A)*c;