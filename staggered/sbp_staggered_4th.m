function [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_4th(n,h,x)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_4th(n,h,x)
% n : Specifies the number of grid points, n+1 (x+ grid) and n+2 (x- grid)
% h : Specifies the grid spacing


assert(n >= 7,'Not enough grid points');  

if nargin < 3 || isempty(x)
  x =[...
     0.054023823103802;  
    -0.937319605810874;  
     1.091344129793494;];
end

pm4  = x(3);
qm03 = x(1);
qm43 = x(2);

% Make sure that the quadratures are positive
pm4 = min(max(0.68552,pm4),1.2171);


% Coefficients determined such that the SBP property is satisfied
qp22 = -227*pm4/21 + 15*qm03/4 + 9*qm43 + 7663/378;
qm31 = 101*pm4/7 - 9*qm03/8 - 9*qm43 - 12487/504;
pp3 = -16*pm4/21 + 2783/1512;
qp02 = -16*pm4/7 + 5*qm03/4 + 3*qm43 + 2951/504;
qm30 = -101*pm4/21 + 3*qm03/8 + 3*qm43 + 3127/378;
qp03 = 101*pm4/21 - 3*qm03/8 - 3*qm43 - 3127/378;
qp31 = -16*pm4/21 + 15*qm03/8 + qm43 + 743/378;
qm21 = -227*pm4/21 + 15*qm03/4 + 9*qm43 + 7663/378;
qp14 = 5*pm4 - 3*qm43 - 67/8;
pm2 = 83*pm4/21 - 4093/1512;
qm10 = 5*pm4/3 + 15*qm03/8 + qm43 - 107/108;
qp30 = -qm03;
pm3 = -23*pm4/7 + 2125/504;
qm20 = 16*pm4/7 - 5*qm03/4 - 3*qm43 - 2951/504;
qm11 = -pm4/7 - 45*qm03/8 - 3*qm43 - 247/84;
pp0 = 16*pm4/21 - 88/189;
qm42 = 3*pm4 - 3*qm43 - 19/3;
qp01 = -5*pm4/3 - 15*qm03/8 - qm43 + 107/108;
pp1 = -16*pm4/7 + 1859/504;
qp12 = 227*pm4/21 - 15*qm03/4 - 9*qm43 - 7663/378;
qm02 = -8*pm4/21 - 3*qm03 + 44/189;
qm22 = 227*pm4/21 - 15*qm03/4 - 9*qm43 - 7663/378;
qp11 = pm4/7 + 45*qm03/8 + 3*qm43 + 247/84;
qm23 = -16*pm4/7 + 5*qm03/4 + 3*qm43 + 2951/504;
qm13 = 16*pm4/21 - 15*qm03/8 - qm43 - 743/378;
qp20 = 8*pm4/21 + 3*qm03 - 44/189;
qm12 = -16*pm4/7 + 45*qm03/8 + 3*qm43 + 743/126;
qp10 = -32*pm4/21 - 3*qm03 + 176/189;
qm40 = 2*pm4 - qm43 - 25/8;
qp04 = -2*pm4 + qm43 + 25/8;
qm32 = -78*pm4/7 + 9*qm03/8 + 9*qm43 + 430/21;
qm41 = -5*pm4 + 3*qm43 + 67/8;
qp21 = 16*pm4/7 - 45*qm03/8 - 3*qm43 - 743/126;
qm33 = 32*pm4/21 - 3*qm03/8 - 3*qm43 - 743/189;
qp33 = -32*pm4/21 + 3*qm03/8 + 3*qm43 + 743/189;
qm01 = 32*pm4/21 + 3*qm03 - 176/189;
qp24 = -3*pm4 + 3*qm43 + 19/3;
qp32 = 16*pm4/7 - 5*qm03/4 - 3*qm43 - 2951/504;
pm1 = -17*pm4/7 + 745/252;
qm00 = -8*pm4/7 - qm03 + 44/63;
qp34 = -qm43;
qp00 = 8*pm4/7 + qm03 - 107/63;
qp13 = -101*pm4/7 + 9*qm03/8 + 9*qm43 + 12487/504;
pm0 = 16*pm4/21 - 88/189;
pp2 = 16*pm4/7 - 197/126;
qp23 = 78*pm4/7 - 9*qm03/8 - 9*qm43 - 430/21;



% Number of coefficients
b = 4;

% Q+ and Q-, top-left corner
QpL = [...
qp00, qp01, qp02, qp03, qp04;
 qp10, qp11, qp12, qp13, qp14;
 qp20, qp21, qp22, qp23, qp24;
 qp30, qp31, qp32, qp33, qp34
];
QmL = [...
qm00, qm01, qm02, qm03;
 qm10, qm11, qm12, qm13;
 qm20, qm21, qm22, qm23;
 qm30, qm31, qm32, qm33;
 qm40, qm41, qm42, qm43
];

% Q+ and Q-
w = b; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qp(end,:) = [];
Qm(:,end) = [];

% Add SBP boundary closures
Qp(1:b,1:b+1) = QpL;
Qp(end-b+1:end,end-b:end) = -fliplr(flipud(QpL));
Qm(1:b+1,1:b) = QmL;
Qm(end-b:end,end-b+1:end) = -fliplr(flipud(QmL));

% P+ and P-
Pp = ones(n+1,1);
Pm = ones(n+2,1);

Pp(1:b) = [pp0,  pp1,  pp2,  pp3]; 
Pp(end-b+1:end) = Pp(b:-1:1);
Pm(1:b+1) = [pm0,  pm1,  pm2,  pm3,  pm4];
Pm(end-b:end) = Pm(b+1:-1:1);
Pp = spdiags(Pp,0,n+1,n+1);
Pm = spdiags(Pm,0,n+2,n+2);

Pp = h*Pp;
Pm = h*Pm;

% nodal and cell-centered grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  
