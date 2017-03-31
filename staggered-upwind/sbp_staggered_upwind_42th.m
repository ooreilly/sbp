function [xp,xm,HIp,HIm,Dp,Dm,Dpp,Dpm,Dmp,Dmm] = sbp_staggered_upwind_42th(n,h)

assert(n >= 8,'Not enough grid points');  

np=n+1;
nm=n+2;


%H- norm
Hm_U=[0.91e2 / 0.720e3 0 0 0; 0 0.325e3 / 0.384e3 0 0; 0 0 0.595e3 / 0.576e3 0; 0 0 0 0.1909e4 / 0.1920e4;];
%H+ norm
Hp_U=[0.805e3 / 0.2304e4 0 0 0; 0 0.955e3 / 0.768e3 0 0; 0 0 0.677e3 / 0.768e3 0; 0 0 0 0.2363e4 / 0.2304e4;];
% Upper part of Qpp. Notice that the B matrix is not included here.
% Qpp+Qpp^T=S/2,Qpp+Qpm^T=0 
Qpp_U=[-0.11e2 / 0.1536e4 0.1493e4 / 0.2304e4 -0.571e3 / 0.4608e4 -0.13e2 / 0.768e3; -0.697e3 / 0.1152e4 -0.121e3 / 0.1536e4 0.97e2 / 0.128e3 -0.473e3 / 0.4608e4; 0.373e3 / 0.4608e4 -0.139e3 / 0.256e3 -0.319e3 / 0.1536e4 0.1999e4 / 0.2304e4; 0.1e1 / 0.32e2 -0.121e3 / 0.4608e4 -0.277e3 / 0.576e3 -0.143e3 / 0.512e3;];

% Upper part of Qmp. Notice that the B matrix is not included here.
% Qmp+Qmp^T=S/2,Qmp+Qmm^T=0 
Qmp_U=[-0.209e3 / 0.5970e4 0.13439e5 / 0.17910e5 -0.13831e5 / 0.47760e5 0.10637e5 / 0.143280e6; -0.44351e5 / 0.71640e5 -0.20999e5 / 0.152832e6 0.230347e6 / 0.229248e6 -0.70547e5 / 0.254720e6; 0.3217e4 / 0.15920e5 -0.86513e5 / 0.114624e6 -0.15125e5 / 0.76416e5 0.1087241e7 / 0.1146240e7; -0.1375e4 / 0.28656e5 0.36117e5 / 0.254720e6 -0.655601e6 / 0.1146240e7 -0.211717e6 / 0.764160e6;];

% The staggered + operator, upper part. Notice that the B matrix is not included here.Qp+Qm^T=0
Qp_U=[-0.338527e6 / 0.1004160e7 0.4197343e7 / 0.4819968e7 0.1423e4 / 0.803328e6 -0.854837e6 / 0.24099840e8; -0.520117e6 / 0.3012480e7 -0.492581e6 / 0.535552e6 0.2476673e7 / 0.2409984e7 0.520117e6 / 0.8033280e7; -0.50999e5 / 0.3012480e7 0.117943e6 / 0.1606656e7 -0.2476673e7 / 0.2409984e7 0.2712193e7 / 0.2677760e7; 0.26819e5 / 0.1004160e7 -0.117943e6 / 0.4819968e7 -0.1423e4 / 0.803328e6 -0.26119411e8 / 0.24099840e8;];

Hp=spdiags(ones(np,1),0,np,np);
Hp(1:4,1:4)=Hp_U;
Hp(np-3:np,np-3:np)=fliplr(flipud(Hp_U));
Hp=Hp*h;
HIp=inv(Hp);

Hm=spdiags(ones(nm,1),0,nm,nm);
Hm(1:4,1:4)=Hm_U;
Hm(nm-3:nm,nm-3:nm)=fliplr(flipud(Hm_U));
Hm=Hm*h;
HIm=inv(Hm);

Qpp=spdiags(repmat([0.7e1 / 0.128e3 -0.67e2 / 0.128e3 -0.55e2 / 0.192e3 0.61e2 / 0.64e2 -0.29e2 / 0.128e3 0.11e2 / 0.384e3;],[np,1]),-2:3,np,np);
Qpp(1:4,1:4)=Qpp_U;
Qpp(np-3:np,np-3:np)=flipud( fliplr(Qpp_U(1:4,1:4) ) )'; 

Qpm=-Qpp';

Qmp=spdiags(repmat([0.7e1 / 0.128e3 -0.67e2 / 0.128e3 -0.55e2 / 0.192e3 0.61e2 / 0.64e2 -0.29e2 / 0.128e3 0.11e2 / 0.384e3;],[nm,1]),-2:3,nm,nm);
Qmp(1:4,1:4)=Qmp_U;
Qmp(nm-3:nm,nm-3:nm)=flipud( fliplr(Qmp_U(1:4,1:4) ) )'; 

Qmm=-Qmp';


Bpp=spalloc(np,np,2);Bpp(1,1)=-1;Bpp(np,np)=1;
Bmp=spalloc(nm,nm,2);Bmp(1,1)=-1;Bmp(nm,nm)=1;

Dpp=HIp*(Qpp+1/2*Bpp) ;
Dpm=HIp*(Qpm+1/2*Bpp) ;


Dmp=HIm*(Qmp+1/2*Bmp) ;
Dmm=HIm*(Qmm+1/2*Bmp) ;


%%% Start with the staggered
Qp=spdiags(repmat([1/24 -9/8 9/8 -1/24],[np,1]),-1:2,np,nm);
Qp(1:4,1:4)=Qp_U;
Qp(np-3:np,nm-3:nm)=flipud( fliplr(-Qp_U(1:4,1:4) ) ); 
Qm=-Qp';

Bp=spalloc(np,nm,2);Bp(1,1)=-1;Bp(np,nm)=1;
Bm=Bp';

Dp=HIp*(Qp+1/2*Bp) ;

Dm=HIm*(Qm+1/2*Bm) ;

% grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  

