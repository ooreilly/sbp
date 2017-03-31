function [xp,xm,HIp,HIm,Dp,Dm,Dpp,Dpm,Dmp,Dmm] = sbp_staggered_upwind_21th(n,h)

assert(n >= 6,'Not enough grid points');  

np=n+1;
nm=n+2;


%H- norm
Hm_U=[0.88829e5 / 0.599299e6 0 0; 0 0.8056187e7 / 0.9588784e7 0; 0 0 0.9700117e7 / 0.9588784e7;];
%H+ norm
Hp_U=[0.7489557e7 / 0.19177568e8 0 0; 0 0.1386089e7 / 0.1198598e7 0; 0 0 0.18276939e8 / 0.19177568e8;];
% Upper part of Qpp. Notice that the B matrix is not included here.
% Qpp+Qpp^T=S/2,Qpp+Qpm^T=0 
Qpp_U=[-0.1e1 / 0.32e2 0.12886609e8 / 0.19177568e8 -0.1349263e7 / 0.9588784e7; -0.10489413e8 / 0.19177568e8 -0.3e1 / 0.16e2 0.16482403e8 / 0.19177568e8; 0.187491e6 / 0.2397196e7 -0.9290815e7 / 0.19177568e8 -0.11e2 / 0.32e2;];
% Upper part of Qmp. Notice that the B matrix is not included here.
% Qmp+Qmp^T=S/2,Qmp+Qmm^T=0 
Qmp_U=[-0.2e1 / 0.21e2 0.12495263e8 / 0.16780372e8 -0.7520839e7 / 0.50341116e8; -0.7700871e7 / 0.16780372e8 -0.31e2 / 0.112e3 0.57771939e8 / 0.67121488e8; 0.2726447e7 / 0.50341116e8 -0.31402783e8 / 0.67121488e8 -0.113e3 / 0.336e3;];
% The staggered + operator, upper part. Notice that the B matrix is not included here.Qp+Qm^T=0
Qp_U=[-0.801195e6 / 0.2397196e7 0.16507959e8 / 0.19177568e8 -0.509615e6 / 0.19177568e8; -0.219745e6 / 0.1198598e7 -0.2112943e7 / 0.2397196e7 0.2552433e7 / 0.2397196e7; 0.42087e5 / 0.2397196e7 0.395585e6 / 0.19177568e8 -0.19909849e8 / 0.19177568e8;];
Hp=spdiags(ones(np,1),0,np,np);
Hp(1:3,1:3)=Hp_U;
Hp(np-2:np,np-2:np)=fliplr(flipud(Hp_U));
Hp=Hp*h;
HIp=inv(Hp);

Hm=spdiags(ones(nm,1),0,nm,nm);
Hm(1:3,1:3)=Hm_U;
Hm(nm-2:nm,nm-2:nm)=fliplr(flipud(Hm_U));
Hm=Hm*h;
HIm=inv(Hm);

Qpp=spdiags(repmat([-0.3e1 / 0.8e1 -0.3e1 / 0.8e1 0.7e1 / 0.8e1 -0.1e1 / 0.8e1;],[np,1]),-1:2,np,np);
Qpp(1:3,1:3)=Qpp_U;
Qpp(np-2:np,np-2:np)=flipud( fliplr(Qpp_U(1:3,1:3) ) )'; 

Qpm=-Qpp';

Qmp=spdiags(repmat([-0.3e1 / 0.8e1 -0.3e1 / 0.8e1 0.7e1 / 0.8e1 -0.1e1 / 0.8e1;],[nm,1]),-1:2,nm,nm);
Qmp(1:3,1:3)=Qmp_U;
Qmp(nm-2:nm,nm-2:nm)=flipud( fliplr(Qmp_U(1:3,1:3) ) )'; 

Qmm=-Qmp';


Bpp=spalloc(np,np,2);Bpp(1,1)=-1;Bpp(np,np)=1;
Bmp=spalloc(nm,nm,2);Bmp(1,1)=-1;Bmp(nm,nm)=1;

Dpp=HIp*(Qpp+1/2*Bpp) ;
Dpm=HIp*(Qpm+1/2*Bpp) ;


Dmp=HIm*(Qmp+1/2*Bmp) ;
Dmm=HIm*(Qmm+1/2*Bmp) ;


%%% Start with the staggered
Qp=spdiags(repmat([-1 1],[np,1]),0:1,np,nm);
Qp(1:3,1:3)=Qp_U;
Qp(np-2:np,nm-2:nm)=flipud( fliplr(-Qp_U(1:3,1:3) ) ); 
Qm=-Qp';

Bp=spalloc(np,nm,2);Bp(1,1)=-1;Bp(np,nm)=1;
Bm=Bp';

Dp=HIp*(Qp+1/2*Bp) ;

Dm=HIm*(Qm+1/2*Bm) ;

% grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  

