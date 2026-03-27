%% finite-difference equation
for n=2:nt-1
    % apply Ricker wave force source
    switch sourcetype
        case 1
            %function s=Ricker(t,fc,shift)
            p(xs,zs,n)=p(xs,zs,n)+rickerp(n-1);
            for i=M+1:nx-M-1
                for j=M+1:nz-M-1
                    sum1=0;sum2=0;
                    for m=1:M/2
                        sum1=sum1+am(m)*(p(i+m,j,n)+p(i-m,j,n)-2*p(i,j,n));
                        sum2=sum2+am(m)*(p(i,j+m,n)+p(i,j-m,n)-2*p(i,j,n));
                    end
                    p(i,j,n+1)=2*p(i,j,n)-p(i,j,n-1)+F(i,j)^2*(sum1+sum2);
                end
            end
        case 2
            uz(xs,zs,n) = uz(xs,zs,n) + rickerf(n-1);   %vertical force.
            for i=M+1:nx-M-1 %need M long grid at boundary for calculating
                for j=M+1:nz-M-1
                    % without applying FD coefficient
                    % l2=F(i,j)^2*(2*ux(i+1,j,n)-2*ux(i,j,n)+ux(i-1,j,n));
                    % l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+1,j+1,n)-uz(i+1,j-1,n)-uz(i-1,j+1,n)+uz(i-1,j-1,n));
                    % l4=F(i,j)^2*gamma(i,j)^2*(ux(i,j+1,n) - 2*ux(i,j,n) + ux(i,j-1,n));
                    % ux(i,j,n+1)=l1+l2+l3+l4;
                    % l1=2*uz(i,j,n) - uz(i,j,n-1) ;
                    % l2=F(i,j)^2*(2*uz(i,j+1,n)-2*uz(i,j,n)+uz(i,j-1,n));
                    % l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+1,j+1,n)-ux(i+1,j-1,n)-ux(i-1,j+1,n)+ux(i-1,j-1,n));
                    % l4=F(i,j)^2*gamma(i,j)^2*(uz(i+1,j,n) - 2*uz(i,j,n) + uz(i-1,j,n));
                    % uz(i,j,n+1)=l1+l2+l3+l4;

                    % with
                    l1x=2*ux(i,j,n) - ux(i,j,n-1) ;
                    l1z=2*uz(i,j,n) - uz(i,j,n-1) ;
                    l2x=0;l2z=0;l4x=0;l4z=0;l3x=0;l3z=0;
                    for m=1:M/2
                        %x
                        l2=F(i,j)^2*(ux(i+m,j,n)-2*ux(i,j,n)+ux(i-m,j,n));
                        l2x=l2x+am(m)*l2;
                        l4=F(i,j)^2*gamma(i,j)^2*(ux(i,j+m,n) - 2*ux(i,j,n) + ux(i,j-m,n));
                        l4x=l4x+am(m)*l4;

                        sl3x=0;
                        for k=1:M/2
                            l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+m,j+k,n)-uz(i-m,j+k,n)-uz(i+m,j-k,n)+uz(i-m,j-k,n));
                            %l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+k,j+m,n)-uz(i-k,j+m,n)-uz(i+k,j-m,n)+uz(i-k,j-m,n));
                            sl3x=sl3x+bm(k)*l3;
                        end
                        l3x=l3x+bm(m)*sl3x;

                        %z
                        l2=F(i,j)^2*(uz(i,j+m,n)-2*uz(i,j,n)+uz(i,j-m,n));
                        l2z=l2z+am(m)*l2;
                        l4=F(i,j)^2*gamma(i,j)^2*(uz(i+m,j,n) - 2*uz(i,j,n) + uz(i-m,j,n));
                        l4z=l4z+am(m)*l4;
                        sl3z=0;
                        for k=1:M/2
                            l3k=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+m,j+k,n)-ux(i+m,j-k,n)-ux(i-m,j+k,n)+ux(i-m,j-k,n));
                            %l3k=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+k,j+m,n)-ux(i+k,j-m,n)-ux(i-k,j+m,n)+ux(i-k,j-m,n));
                            sl3z=sl3z+bm(k)*l3k;
                        end
                        l3z=l3z+bm(m)*sl3z;
                    end
                    ux(i,j,n+1)=l1x+l2x+l3x+l4x;
                    uz(i,j,n+1)=l1z+l2z+l3z+l4z;
                end
            end

        case 3
            ux(xs,zs,n)=ux(xs,zs,n)+rickerf(n-1);    %Horizontal force.
            for i=M+1:nx-M-1 %need M long grid at boundary for calculating
                for j=M+1:nz-M-1
                    % without applying FD coefficient
                    % l2=F(i,j)^2*(2*ux(i+1,j,n)-2*ux(i,j,n)+ux(i-1,j,n));
                    % l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+1,j+1,n)-uz(i+1,j-1,n)-uz(i-1,j+1,n)+uz(i-1,j-1,n));
                    % l4=F(i,j)^2*gamma(i,j)^2*(ux(i,j+1,n) - 2*ux(i,j,n) + ux(i,j-1,n));
                    % ux(i,j,n+1)=l1+l2+l3+l4;
                    % l1=2*uz(i,j,n) - uz(i,j,n-1) ;
                    % l2=F(i,j)^2*(2*uz(i,j+1,n)-2*uz(i,j,n)+uz(i,j-1,n));
                    % l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+1,j+1,n)-ux(i+1,j-1,n)-ux(i-1,j+1,n)+ux(i-1,j-1,n));
                    % l4=F(i,j)^2*gamma(i,j)^2*(uz(i+1,j,n) - 2*uz(i,j,n) + uz(i-1,j,n));
                    % uz(i,j,n+1)=l1+l2+l3+l4;

                    % with
                    l1x=2*ux(i,j,n) - ux(i,j,n-1) ;
                    l1z=2*uz(i,j,n) - uz(i,j,n-1) ;
                    l2x=0;l2z=0;l4x=0;l4z=0;l3x=0;l3z=0;
                    for m=1:M/2
                        %x
                        l2=F(i,j)^2*(ux(i+m,j,n)-2*ux(i,j,n)+ux(i-m,j,n));
                        l2x=l2x+am(m)*l2;
                        l4=F(i,j)^2*gamma(i,j)^2*(ux(i,j+m,n) - 2*ux(i,j,n) + ux(i,j-m,n));
                        l4x=l4x+am(m)*l4;
                        %Here also can use 1 m loop for faster calculation.
                        sl3x=0;
                        for k=1:M/2
                            l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+m,j+k,n)-uz(i-m,j+k,n)-uz(i+m,j-k,n)+uz(i-m,j-k,n));
                            %l3=F(i,j)^2/4*(1-gamma(i,j)^2)*(uz(i+k,j+m,n)-uz(i-k,j+m,n)-uz(i+k,j-m,n)+uz(i-k,j-m,n));
                            sl3x=sl3x+bm(k)*l3;
                        end
                        l3x=l3x+bm(m)*sl3x;

                        %z
                        l2=F(i,j)^2*(uz(i,j+m,n)-2*uz(i,j,n)+uz(i,j-m,n));
                        l2z=l2z+am(m)*l2;
                        l4=F(i,j)^2*gamma(i,j)^2*(uz(i+m,j,n) - 2*uz(i,j,n) + uz(i-m,j,n));
                        l4z=l4z+am(m)*l4;
                        sl3z=0;
                        for k=1:M/2
                            l3k=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+m,j+k,n)-ux(i+m,j-k,n)-ux(i-m,j+k,n)+ux(i-m,j-k,n));
                            %l3k=F(i,j)^2/4*(1-gamma(i,j)^2)*(ux(i+k,j+m,n)-ux(i+k,j-m,n)-ux(i-k,j+m,n)+ux(i-k,j-m,n));
                            sl3z=sl3z+bm(k)*l3k;
                        end
                        l3z=l3z+bm(m)*sl3z;
                    end
                    ux(i,j,n+1)=l1x+l2x+l3x+l4x;
                    uz(i,j,n+1)=l1z+l2z+l3z+l4z;
                end
            end
    end
end