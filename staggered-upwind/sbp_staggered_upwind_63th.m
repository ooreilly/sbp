function [xp,xm,HIp,HIm,Dp,Dm,Dpp,Dpm,Dmp,Dmm] = sbp_staggered_upwind_63th(n,h)

assert(n >= 12,'Not enough grid points');  

np=n+1;
nm=n+2;

bp=6; % Number of boundary points


%H- norm
Hm_U=[0.18373e5 / 0.136080e6 0 0 0 0 0; 0 0.228247e6 / 0.276480e6 0 0 0 0; 0 0 0.219557e6 / 0.207360e6 0 0 0; 0 0 0 0.44867e5 / 0.46080e5 0 0; 0 0 0 0 0.487757e6 / 0.483840e6 0; 0 0 0 0 0 0.2485447e7 / 0.2488320e7;];
%H+ norm
Hp_U=[0.174529e6 / 0.552960e6 0 0 0 0 0; 0 0.769723e6 / 0.552960e6 0 0 0 0; 0 0 0.172613e6 / 0.276480e6 0 0 0; 0 0 0 0.343867e6 / 0.276480e6 0 0; 0 0 0 0 0.503237e6 / 0.552960e6 0; 0 0 0 0 0 0.560831e6 / 0.552960e6;];
% Upper part of Qpp. Notice that the B matrix is not included here.
% Qpp+Qpp^T=S/2,Qpp+Qpm^T=0 
Qpp_U=[-0.11e2 / 0.12288e5 0.3248629e7 / 0.4976640e7 -0.742121e6 / 0.9953280e7 -0.173089e6 / 0.1658880e7 0.100627e6 / 0.9953280e7 0.263e3 / 0.15552e5; -0.3212989e7 / 0.4976640e7 -0.341e3 / 0.20480e5 0.488429e6 / 0.995328e6 0.2108717e7 / 0.9953280e7 0.5623e4 / 0.829440e6 -0.468779e6 / 0.9953280e7; 0.635201e6 / 0.9953280e7 -0.2135641e7 / 0.4976640e7 -0.2233e4 / 0.30720e5 0.298757e6 / 0.497664e6 -0.2148901e7 / 0.9953280e7 0.99581e5 / 0.1658880e7; 0.184969e6 / 0.1658880e7 -0.2671829e7 / 0.9953280e7 -0.1044721e7 / 0.2488320e7 -0.4697e4 / 0.30720e5 0.882599e6 / 0.995328e6 -0.2114063e7 / 0.9953280e7; -0.118447e6 / 0.9953280e7 0.15761e5 / 0.829440e6 0.915757e6 / 0.9953280e7 -0.2923243e7 / 0.4976640e7 -0.4279e4 / 0.20480e5 0.4614511e7 / 0.4976640e7; -0.263e3 / 0.15552e5 0.422447e6 / 0.9953280e7 -0.5185e4 / 0.331776e6 0.424727e6 / 0.9953280e7 -0.570779e6 / 0.995328e6 -0.2761e4 / 0.12288e5;];

% Upper part of Qmp. Notice that the B matrix is not included here.
% Qmp+Qmp^T=S/2,Qmp+Qmm^T=0 
Qmp_U=[-0.6660404399e10 / 0.975680535935e12 0.42802131970831759e17 / 0.60695134779444480e17 -0.27241603626152813e17 / 0.91042702169166720e17 0.148772145985039e15 / 0.1264481974571760e16 -0.386865438537449e15 / 0.30347567389722240e17 -0.739491571084877e15 / 0.182085404338333440e18; -0.40937739835891039e17 / 0.60695134779444480e17 -0.20934998251893e14 / 0.570912496455680e15 0.3069347264824655e16 / 0.3082927480860672e16 -0.531249064089227e15 / 0.1541463740430336e16 0.6343721417240047e16 / 0.107902461830123520e18 0.194756943144697e15 / 0.138731736638730240e18; 0.24212381292051773e17 / 0.91042702169166720e17 -0.13932096926615587e17 / 0.15414637404303360e17 -0.22149140293797e14 / 0.285456248227840e15 0.443549615179363e15 / 0.481707418884480e15 -0.1531771257444041e16 / 0.5994581212784640e16 0.23579222779798361e17 / 0.416195209916190720e18; -0.119651391006031e15 / 0.1264481974571760e16 0.2061363806050549e16 / 0.7707318702151680e16 -0.355212104634871e15 / 0.481707418884480e15 -0.42909037900311e14 / 0.285456248227840e15 0.25466291778687943e17 / 0.26975615457530880e17 -0.3289076301679109e16 / 0.11560978053227520e17; 0.153994229603129e15 / 0.30347567389722240e17 -0.2655886084488631e16 / 0.107902461830123520e18 0.782300684927837e15 / 0.5994581212784640e16 -0.17511430871269903e17 / 0.26975615457530880e17 -0.413193098349471e15 / 0.1998193737594880e16 0.189367309285289755e18 / 0.194224431294222336e18; 0.894580992211517e15 / 0.182085404338333440e18 -0.209441083772219e15 / 0.27746347327746048e17 -0.4946149632449393e16 / 0.416195209916190720e18 0.334964642443661e15 / 0.2890244513306880e16 -0.604469352802317407e18 / 0.971122156471111680e18 -0.128172128502407e15 / 0.570912496455680e15;];

% The staggered + operator, upper part. Notice that the B matrix is not included here.Qp+Qm^T=0
Qp_U=[-0.34660470729017653729e20 / 0.113641961250214656000e21 0.351671379135966469961e21 / 0.415604886857927884800e21 0.1819680091728191503e19 / 0.103901221714481971200e21 -0.18252344147469061739e20 / 0.346337405714939904000e21 -0.18145368485798816351e20 / 0.727308552001373798400e21 0.2627410615589536403e19 / 0.138534962285975961600e21; -0.1606450873889019037e19 / 0.7576130750014310400e19 -0.8503979509850519441e19 / 0.9235664152398397440e19 0.2208731907526094393e19 / 0.2308916038099599360e19 0.1143962309827873891e19 / 0.7696386793665331200e19 0.1263616990270014071e19 / 0.16162412266697195520e20 -0.1402288892096389187e19 / 0.27706992457195192320e20; -0.502728075208147729e18 / 0.11364196125021465600e20 0.5831273443201206481e19 / 0.41560488685792788480e20 -0.9031420599281409001e19 / 0.10390122171448197120e20 0.29977986617775158621e20 / 0.34633740571493990400e20 -0.7995649008389734663e19 / 0.72730855200137379840e20 0.728315692435313537e18 / 0.41560488685792788480e20; 0.710308100786581369e18 / 0.11364196125021465600e20 -0.2317346723533341809e19 / 0.41560488685792788480e20 -0.1357359577229545879e19 / 0.10390122171448197120e20 -0.35124499190079631261e20 / 0.34633740571493990400e20 0.81675241511291974823e20 / 0.72730855200137379840e20 0.144034831596315317e18 / 0.13853496228597596160e20; 0.13360631165154733e17 / 0.841792305557145600e18 -0.875389186128426797e18 / 0.27706992457195192320e20 0.95493318392786453e17 / 0.6926748114298798080e19 0.1714625642820850967e19 / 0.23089160380995993600e20 -0.53371483072841696197e20 / 0.48487236800091586560e20 0.30168063964639488547e20 / 0.27706992457195192320e20; -0.1943232250834614071e19 / 0.113641961250214656000e21 0.8999269402554660119e19 / 0.415604886857927884800e21 0.1242786058815312817e19 / 0.103901221714481971200e21 -0.7480218714053301461e19 / 0.346337405714939904000e21 0.28468183824757664191e20 / 0.727308552001373798400e21 -0.476082521721490837529e21 / 0.415604886857927884800e21;];

Hp=spdiags(ones(np,1),0,np,np);
Hp(1:bp,1:bp)=Hp_U;
Hp(np-bp+1:np,np-bp+1:np)=fliplr(flipud(Hp_U));
Hp=Hp*h;
HIp=inv(Hp);

Hm=spdiags(ones(nm,1),0,nm,nm);
Hm(1:bp,1:bp)=Hm_U;
Hm(nm-bp+1:nm,nm-bp+1:nm)=fliplr(flipud(Hm_U));
Hm=Hm*h;
HIm=inv(Hm);

Qpp=spdiags(repmat([-0.157e3 / 0.15360e5 0.537e3 / 0.5120e4 -0.3147e4 / 0.5120e4 -0.231e3 / 0.1024e4 0.999e3 / 0.1024e4 -0.1461e4 / 0.5120e4 0.949e3 / 0.15360e5 -0.33e2 / 0.5120e4;],[np,1]),-3:4,np,np);
Qpp(1:bp,1:bp)=Qpp_U;
Qpp(np-bp+1:np,np-bp+1:np)=flipud( fliplr(Qpp_U ) )'; 

Qpm=-Qpp';

Qmp=spdiags(repmat([-0.157e3 / 0.15360e5 0.537e3 / 0.5120e4 -0.3147e4 / 0.5120e4 -0.231e3 / 0.1024e4 0.999e3 / 0.1024e4 -0.1461e4 / 0.5120e4 0.949e3 / 0.15360e5 -0.33e2 / 0.5120e4;],[nm,1]),-3:4,nm,nm);
Qmp(1:bp,1:bp)=Qmp_U;
Qmp(nm-bp+1:nm,nm-bp+1:nm)=flipud( fliplr(Qmp_U ) )'; 

Qmm=-Qmp';

Bpp=spalloc(np,np,2);Bpp(1,1)=-1;Bpp(np,np)=1;
Bmp=spalloc(nm,nm,2);Bmp(1,1)=-1;Bmp(nm,nm)=1;


Dpp=HIp*(Qpp+1/2*Bpp) ;
Dpm=HIp*(Qpm+1/2*Bpp) ;


Dmp=HIm*(Qmp+1/2*Bmp) ;
Dmm=HIm*(Qmm+1/2*Bmp) ;


%%% Start with the staggered
Qp=spdiags(repmat([-0.3e1 / 0.640e3 0.25e2 / 0.384e3 -0.75e2 / 0.64e2 0.75e2 / 0.64e2 -0.25e2 / 0.384e3 0.3e1 / 0.640e3],[np,1]),-2:3,np,nm);
Qp(1:bp,1:bp)=Qp_U;
Qp(np-bp+1:np,nm-bp+1:nm)=flipud( fliplr(-Qp_U ) ); 
Qm=-Qp';

Bp=spalloc(np,nm,2);Bp(1,1)=-1;Bp(np,nm)=1;
Bm=Bp';

Dp=HIp*(Qp+1/2*Bp) ;

Dm=HIm*(Qm+1/2*Bm) ;

% grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  
