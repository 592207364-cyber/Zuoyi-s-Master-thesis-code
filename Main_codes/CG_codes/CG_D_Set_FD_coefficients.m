%% get am  a0 
m2k=zeros(M/2);
for n=1:M/2
    for i=1:M/2
        m2k(n,i)= i^(2*(n));%i->m n->k
    end
end
Right=zeros(M/2,1);
Right(1,1)=1;
am=m2k\Right;
% a0=2*sum(am);

%% get bm
m2k_1=zeros(M/2);
for n=1:M/2
    for i=1:M/2
        m2k_1(n,i)=i^(2*(n)-1);
    end
end
bm=m2k_1\Right;