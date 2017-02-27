function [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered(order, n, h)
% [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered(order, n, h)
% Construct SBP staggered grid operators
%
% Input arguments:
% order : Order of accuracy in the interior (2, 4, or 6)
% n     : Number of grid points n+1 (x+ grid) n+2 (x- grid)
% h     : Grid spacing
%
% Output arguments:
% xp, xm           : Grids 
% Pp, Pm, Qp, Qm,  : SBP staggered grid operators

    switch order
      case 2
        [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered_2(n, h);
      case 4
        [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered_4(n, h);
      case 6
        [xp, xm, Pp, Pm, Qp, Qm] = sbp_staggered_6(n, h);
      otherwise
       error('SBP staggered grid operator not implemented');   
      end

xp = xp';
xm = xm';

end
