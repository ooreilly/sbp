function [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered(order,n,h,truncate,test)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered(order,n,h,truncate,test)
%
% Input:
% order                : Order of accuracy (2,4,6, or 8)
% n                    : Number of grid points xp (n+1) 
% h                    : Grid spacing
% truncate (true/false): Truncate first row and last row of Dm
%
% Optional:
%  test (true/false) : Tests that the SBP operators are accurate.
%
% Output:
% xp                 : Grid defined by x_j = h*j
% xm                 : Grid defined by x_j = h*(j+1/2) in the interior
% Pp, Pm             : Quadrature rules
% Qp, Qm             : Operators that satisfy SBP property Qp + Qm^T = B

 if nargin < 5
   test = false;
 end

 switch order
   case 2
     [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_2nd(n,h,test);
   case 4
     [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_4th(n,h,test);
   case 6
     [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_6th(n,h,test);
   case 8
     [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_8th(n,h,test);
   otherwise
    error('SBP staggered grid injection operator not implemented');   
 end

 if truncate
   Pm(:,1) = [];
   Pm(:,end) = [];
   Pm(1,:) = [];
   Pm(end,:) = [];
   Qp(:,1) = [];
   Qp(:,end) = [];
   Qm(1,:) = [];
   Qm(end,:) = [];   
   xm = xm(2:end-1);
 end

end
