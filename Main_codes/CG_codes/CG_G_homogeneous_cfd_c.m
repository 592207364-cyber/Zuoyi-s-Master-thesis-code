kh=0.01:0.01:pi;
if sourcetype==1
    [cgcfdvp_cacoustic,~]=Z_Function_CG_homogeneous_cfd_c(phi,vp,dt,dh,M,kh,sourcetype,am,bm,vs0/vp0);
else 
    [cgcfdvp_c,cgcfdvs_c]=Z_Function_CG_homogeneous_cfd_c(phi,vp,dt,dh,M,kh,sourcetype,am,bm,gamma0);
end