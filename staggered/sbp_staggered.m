function [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered(order,n,h,x)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered(order,n,h,x)
% Construct SBP staggered grid operators
%
% Input:
% order : Specifies the order of accuracy in the interior
% n     : Specifies the number of grid points n+1 (x+ grid) n+2 (x- grid)
% h     : Specifies the grid spacing
%
% Optional:
% x    : A vector of length m, specifying the free parameters in the boundary closure. 
%        These parameters can be selected to minimize for example:
%        the truncation error or spectral radius of the operators. Change it if you want to override
%        the optimized values provided. 
%        Number of free parameters:
%        order   m
%        2       0
%        4       3
%        6       7
%
% Output
% xp,xm         : Grids
% Pp,Pm,Qp,Qm,  : SBP staggered grid operators

if nargin < 4
  test = false;
end
if nargin < 5
  x = [];
end

    switch order
      case 2
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_2nd(n,h);
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_4th(n,h,x);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_6th(n,h,x);
      otherwise
       error('SBP staggered grid operator not implemented');   
      end

xp = xp';
xm = xm';

end
