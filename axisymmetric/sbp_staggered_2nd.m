function [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_2nd(n,h,test)

if nargin < 3
  test = false;
end

assert(n >= 4,'Not enough grid points');

% Coefficients determined such that the SBP property is satisfied
pm1 = 1;
qm01 = -pm1 + 2;
qp11 = pm1 - 2;
qm11 = pm1 - 2;
qp02 = pm1 - 1;
pp0 = pm1 - 1/2;
qm00 = pm1 - 2;
qp10 = 0;
qp12 = -pm1 + 2;
qp00 = -1;
pp1 = -pm1 + 2;
qp01 = -pm1 + 2;
pm0 = -pm1 + 2;
qm10 = -pm1 + 1;



% Number of coefficients
b = 2;

% Q+ and Q-, top-left corner
QpL = [...
qp00, qp01, qp02;
 qp10, qp11, qp12
];
QmL = [...
0, 0;
 qm00, qm01;
 qm10, qm11
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

Pp(1:b) = [pp0,  pp1]; 
Pp(end-b+1:end) = Pp(b:-1:1);
Pm(1:b+1) = [0,  pm0,  pm1];
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
