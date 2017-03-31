function [x, xhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat, HI, HIhat] = sbp_staggered_upwind(p, n, h)
% [x, xhat, HI, HIhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat] = sbp_staggered_upwind(p, n, h)
%
% This function constructs SBP staggered and upwind operators of the first derivative.
% The order of accuracy of the operators is 'p/2' on the boundary and 'p' in the interior.
% All of the operators are constructed using two different grids, a regular grid 'x' and a non-regular grid 'xhat'.
% The 'hat' always refers to the non-regular grid.
% To weakly impose boundary conditions using SAT terms, the inverse norms 'HI' and 'HIhat' are also provided.  
%
% Input arguments:
%   p     : Specifies the order of accuracy in the interior.
%   n     : Specifies the number of grid points.
%   h     : Specifies the grid spacing.
%
% Output arguments:
% x             : Regular grid,         x = [0, h,  ...,  n*h]            (n+1 grid points).
% xhat          : Non-regular grid,  xhat = [0,h/2, ..., (n-1/2)*h, n*h]  (n+2 grid points).
% Ds, Dshat     : Staggered grid operators.
% Dp, Dm        : Upwind operators, regular grid.
% Dphat, Dmhat  : Upwind operators, non-regular grid.
% HI, HIhat     : Inverse norms.

    switch p
      case 2
        [x, xhat, HI, HIhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat] = sbp_staggered_upwind_21th(n,h);
      case 4
        [x, xhat, HI, HIhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat] = sbp_staggered_upwind_42th(n,h);
      case 6
        [x, xhat, HI, HIhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat] = sbp_staggered_upwind_63th(n,h);
      case 8
        [x, xhat, HI, HIhat, Ds, Dshat, Dp, Dm, Dphat, Dmhat] = sbp_staggered_upwind_84th(n,h);
      otherwise
       error('SBP staggered-upwind grid operator not implemented');   
      end

x    = x';
xhat = xhat';


end
