function [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_4nd(n,h,test)

if nargin < 3
  test = false;
end

assert(n >= 8,'Not enough grid points');

% Free parameters determined by optimizing spectral radius/truncation error
x = [0.002435378086542   0.999421296296229];
qm03 = x(1);
pm3 = x(2);


% Coefficients determined such that the SBP property is satisfied
qp22 = -3*pm3 - 9*qm03 + 2;
qm31 = -2*pm3 - 3*qm03 + 49/24;
pp3 = -pm3 + 143/72;
qp02 = -3*qm03 - 1/24;
qm30 = pm3 + qm03 - 73/72;
qp03 = 2*pm3 + 3*qm03 - 37/18;
qp31 = -qm03;
qm21 = 6*pm3 + 9*qm03 - 49/8;
qp14 = 2*pm3 + 3*qm03 - 49/24;
pm2 = -3*pm3 + 97/24;
qm10 = 3*qm03 + 1/24;
qp30 = 0;
qm20 = -2*pm3 - 3*qm03 + 37/18;
qm11 = -3*pm3 - 9*qm03 + 2;
pp0 = pm3 - 11/18;
qp01 = -pm3 + qm03 + 25/12;
pp1 = -3*pm3 + 33/8;
qp12 = 3*pm3 + 9*qm03 - 2;
qm02 = -3*qm03;
qm22 = -3*pm3 - 9*qm03 + 2;
qp11 = pm3 - 3*qm03 - 25/12;
qm23 = -pm3 + 3*qm03 + 19/9;
qm13 = -3*qm03 - 1/24;
qp20 = 0;
qm12 = 3*pm3 + 9*qm03 - 2;
qp10 = 0;
qp04 = -pm3 - qm03 + 73/72;
qm32 = 3*qm03;
qp21 = 3*qm03;
qm33 = pm3 - qm03 - 19/9;
qp33 = pm3 - 3*qm03 - 19/9;
qm01 = -pm3 + 3*qm03 + 25/12;
qp24 = -3*qm03;
qp32 = 3*qm03 + 1/24;
pm1 = 3*pm3 - 17/8;
qm00 = pm3 - qm03 - 25/12;
qp34 = -pm3 + qm03 + 19/9;
qp00 = -1;
qp13 = -6*pm3 - 9*qm03 + 49/8;
pm0 = -pm3 + 25/12;
pp2 = 3*pm3 - 2;
qp23 = 3*pm3 + 9*qm03 - 2;



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
0, 0, 0, 0;
 qm00, qm01, qm02, qm03;
 qm10, qm11, qm12, qm13;
 qm20, qm21, qm22, qm23;
 qm30, qm31, qm32, qm33
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
Pm(1:b+1) = [0,  pm0,  pm1,  pm2,  pm3];
Pm(end-b:end) = Pm(b+1:-1:1);
Pp = spdiags(Pp,0,n+1,n+1);
Pm = spdiags(Pm,0,n+2,n+2);

Pp = h*Pp;
Pm = h*Pm;

xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  


% Test operators
if test
for j=0:b/2
  disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
  num2str(max(abs(Qp*xm.^j-j*Pp*xp.^max([j-1,0]))))]);
  disp([ 'Dm, j = ' num2str(j) ' Error max = '...
  num2str(max(abs(Qm*xp.^j-j*Pm*xm.^max([j-1,0]))))]);
end  
end
