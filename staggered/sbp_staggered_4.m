function [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered_4(n, h)
% [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered_4(n, h)
% Input arguments:
% n : number of grid points.   (n+1) xp grid (n+2) xm grid
% h : grid spacing


assert(n >= 7, 'Not enough grid points');  

% Q+ and Q-, top-left corner
QpL = [...
-323/378, 2783/3072, -1085/27648, -17/9216, -667/64512;
 -5/21, -847/1024, 1085/1024, -37/1024, 899/21504;
 5/42, -121/1024, -1085/1024, 3575/3072, -753/7168;
 -5/189, 121/3072, 1085/27648, -10759/9216, 74635/64512
];
QmL = [...
-55/378, 5/21, -5/42, 5/189;
 -2783/3072, 847/1024, 121/1024, -121/3072;
 1085/27648, -1085/1024, 1085/1024, -1085/27648;
 17/9216, 37/1024, -3575/3072, 10759/9216;
 667/64512, -899/21504, 753/7168, -74635/64512
];

% Q+ and Q-
w = 4; 
s = [ 1/24,  -9/8,  9/8,  -1/24];  
Qp = spdiags(repmat(-s(end:-1:1),[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qp(end,:) = [];
Qm(:,end) = [];

% Add SBP boundary closures
bp = 4; 
bm = 5;
Qp(1:bp,1:bm) = double(QpL);
Qp(end-bp+1:end,end-bm+1:end) = -fliplr(flipud(QpL));
Qm(1:bm,1:bp) = double(QmL);
Qm(end-bm+1:end,end-bp+1:end) = -fliplr(flipud(QmL));

% P+ and P-
Pp = ones(n+1,1);
Pm = ones(n+2,1);

Pp(1:bp) = [407/1152,  473/384,  343/384,  1177/1152]; 
Pp(end-bp+1:end) = Pp(bp:-1:1);
Pm(1:bm) = [5/63,  121/128,  1085/1152,  401/384,  2659/2688];
Pm(end-bm+1:end) = Pm(bm:-1:1);
Pp = spdiags(Pp,0,n+1,n+1);
Pm = spdiags(Pm,0,n+2,n+2);

Pp = h*Pp;
Pm = h*Pm;

xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  
