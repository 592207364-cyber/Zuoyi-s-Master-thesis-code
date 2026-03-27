% Loop over time steps and grid points
% time loop in Matlab is dt to t.
switch sourcetype
    case 1
        for n=1:nt-1
            sxx(xs,zs,n)=sxx(xs,zs,n)+rickerp(n);%%function s=Ricker(t,fc,shift)
            szz(xs,zs,n)=szz(xs,zs,n)+rickerp(n);%% s= source weight* source wavelet
            % in ppt source should add at sxx n^+  (1.5,...) so here shift +0.5
            %consider first loop: code_sxx 1 ->sxx 0.5
            %when we have sxx n+(0.5) with  vx vz n (1) we can get sxx n++1(1.5)
            %and  from vx vz n(1) and sxx szz n++1(1.5),vx vz n+1(2) can be get
            for i=M/2+1:nx-M/2-1
                for j=M/2+1:nz-M/2-1
                    vxx=0.0;vxz=0.0; vzx=0.0;vzz=0.0;
                    for m=1:M/2
                        vxx=vxx+b_m(m)*(vx(i+m-1,j,n)-vx(i-m,j,n));
                        vzx=vzx+b_m(m)*(vz(i+m,j,n)-vz(i-m+1,j,n));
                        vxz=vxz+b_m(m)*(vx(i,j+m,n)-vx(i,j-m+1,n));
                        vzz=vzz+b_m(m)*(vz(i,j+m-1,n)-vz(i,j-m,n));
                    end
                    sxx(i,j,n+1)=sxx(i,j,n)+lambda2mu(i,j)*dt/dh*vxx+lambdadt_dh(i,j)*vzz;%EQ24-1
                    szz(i,j,n+1)=szz(i,j,n)+lambda2mu(i,j)*dt/dh*vzz+lambdadt_dh(i,j)*vxx;%EQ24-2
                    sxz(i,j,n+1)=sxz(i,j,n)+mu(i,j)*dt/dh*(vzx+vxz);%EQ24-3
                    %consider first loop: sxx 2-> sxx 1,5 = sxx 1(0,5)+ vx,vz 1
                    p(i,j,n)=sxx(i,j,n)+szz(i,j,n);
                    %now we have sxx szz sxz 1.5 we need vx vz 1 to get vx vz 2
                end
            end
            %equation of motion
            for i=M/2+1:nx-M/2-1
                for j=M/2+1:nz-M/2-1
                    sxx_x=0;sxz_x=0;sxz_z=0.0;szz_z=0.0;
                    for m=1:M/2
                        %eq25-1
                        sxx_x=sxx_x+b_m(m)*(sxx(i+m,j,n+1)-sxx(i-m+1,j,n+1));
                        sxz_x=sxz_x+b_m(m)*(sxz(i+m-1,j,n+1)-sxz(i-m,j,n+1));
                        %eq25-2
                        sxz_z=sxz_z+b_m(m)*(sxz(i,j+m-1,n+1)-sxz(i,j-m,n+1));
                        szz_z=szz_z+b_m(m)*(szz(i,j+m,n+1)-szz(i,j-m+1,n+1));
                    end
                    vz(i,j,n+1)=vz(i,j,n)+dt_rhoh(i,j)*(sxz_x+szz_z); %ppt L4P21EQ25-1
                    vx(i,j,n+1)=vx(i,j,n)+dt_rhoh(i,j)*(sxx_x+sxz_z); %ppt L4P21EQ25-2
                    %consider first loop: vx 2= vx 1 + sxx szz 1.5
                end
            end
        end
    case 2
        for n=1:nt-1
            vz(xs,zs,n)=vz(xs,zs,n)+rickerf(n);    %vertical force.
            % in ppt source should add at vx n so here no shift
            % first loop: we have vz 1   need sxx syy sxy 0.5 (n,1)
                for i=M/2+1:nx-M/2-1
                    for j=M/2+1:nz-M/2-1
                        vxx=0.0;vxz=0.0; vzx=0.0;vzz=0.0;
                        for m=1:M/2
                            vxx=vxx+b_m(m)*(vx(i+m-1,j,n)-vx(i-m,j,n));
                            vzx=vzx+b_m(m)*(vz(i+m,j,n)-vz(i-m+1,j,n));
                            vxz=vxz+b_m(m)*(vx(i,j+m,n)-vx(i,j-m+1,n));
                            vzz=vzz+b_m(m)*(vz(i,j+m-1,n)-vz(i,j-m,n));
                        end
                        sxx(i,j,n+1)=sxx(i,j,n)+lambda2mu(i,j)*dt/dh*vxx+lambdadt_dh(i,j)*vzz;%EQ24-1
                        szz(i,j,n+1)=szz(i,j,n)+lambda2mu(i,j)*dt/dh*vzz+lambdadt_dh(i,j)*vxx;%EQ24-2
                        sxz(i,j,n+1)=sxz(i,j,n)+mu(i,j)*dt/dh*(vzx+vxz);%EQ24-3
                        %consider first loop: sxx 2-> sxx 1,5 = sxx 1(0,5)+ vx,vz 1
                         p(i,j,n+2)=sxx(i,j,n+1)+szz(i,j,n+1);

                    end
                end
                %equation of motion
                for i=M/2+1:nx-M/2-1
                    for j=M/2+1:nz-M/2-1
                        sxx_x=0;sxz_x=0;sxz_z=0.0;szz_z=0.0;
                        for m=1:M/2
                            %eq25-1
                            sxx_x=sxx_x+b_m(m)*(sxx(i+m,j,n+1)-sxx(i-m+1,j,n+1));
                            sxz_x=sxz_x+b_m(m)*(sxz(i+m-1,j,n+1)-sxz(i-m,j,n+1));
                            %eq25-2
                            sxz_z=sxz_z+b_m(m)*(sxz(i,j+m-1,n+1)-sxz(i,j-m,n+1));
                            szz_z=szz_z+b_m(m)*(szz(i,j+m,n+1)-szz(i,j-m+1,n+1));
                        end
                        vz(i,j,n+1)=vz(i,j,n)+dt_rhoh(i,j)*(sxz_x+szz_z); %ppt L4P21EQ25-1
                        vx(i,j,n+1)=vx(i,j,n)+dt_rhoh(i,j)*(sxx_x+sxz_z); %ppt L4P21EQ25-2
                        %consider first loop: vx 2= vx 1 + sxx szz 1.5
                    end
                end

            end
        case 3
            for n=1:nt-1
                vx(xs,zs,n)=vx(xs,zs,n)+rickerf(n);    %Horizontal force.
            % in ppt source should add at vx n so here no shift
            % first loop: we have vz 1   need sxx syy sxy 0.5 (n,1)
                for i=M/2+1:nx-M/2-1
                    for j=M/2+1:nz-M/2-1
                        vxx=0.0;vxz=0.0; vzx=0.0;vzz=0.0;
                        for m=1:M/2
                            vxx=vxx+b_m(m)*(vx(i+m-1,j,n)-vx(i-m,j,n));
                            vzx=vzx+b_m(m)*(vz(i+m,j,n)-vz(i-m+1,j,n));
                            vxz=vxz+b_m(m)*(vx(i,j+m,n)-vx(i,j-m+1,n));
                            vzz=vzz+b_m(m)*(vz(i,j+m-1,n)-vz(i,j-m,n));
                        end
                        sxx(i,j,n+1)=sxx(i,j,n)+lambda2mu(i,j)*dt/dh*vxx+lambdadt_dh(i,j)*vzz;%EQ24-1
                        szz(i,j,n+1)=szz(i,j,n)+lambda2mu(i,j)*dt/dh*vzz+lambdadt_dh(i,j)*vxx;%EQ24-2
                        sxz(i,j,n+1)=sxz(i,j,n)+mu(i,j)*dt/dh*(vzx+vxz);%EQ24-3
                        %consider first loop: sxx 2-> sxx 1,5 = sxx 1(0,5)+ vx,vz 1
                         p(i,j,n+2)=sxx(i,j,n+1)+szz(i,j,n+1);

                    end
                end
                %equation of motion
                for i=M/2+1:nx-M/2-1
                    for j=M/2+1:nz-M/2-1
                        sxx_x=0;sxz_x=0;sxz_z=0.0;szz_z=0.0;
                        for m=1:M/2
                            %eq25-1
                            sxx_x=sxx_x+b_m(m)*(sxx(i+m,j,n+1)-sxx(i-m+1,j,n+1));
                            sxz_x=sxz_x+b_m(m)*(sxz(i+m-1,j,n+1)-sxz(i-m,j,n+1));
                            %eq25-2
                            sxz_z=sxz_z+b_m(m)*(sxz(i,j+m-1,n+1)-sxz(i,j-m,n+1));
                            szz_z=szz_z+b_m(m)*(szz(i,j+m,n+1)-szz(i,j-m+1,n+1));
                        end
                        vz(i,j,n+1)=vz(i,j,n)+dt_rhoh(i,j)*(sxz_x+szz_z); %ppt L4P21EQ25-1
                        vx(i,j,n+1)=vx(i,j,n)+dt_rhoh(i,j)*(sxx_x+sxz_z); %ppt L4P21EQ25-2
                        %consider first loop: vx 2= vx 1 + sxx szz 1.5
                    end
                end
            end
end