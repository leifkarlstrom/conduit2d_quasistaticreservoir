function [xp,xm,Pp,Pm,Qp,Qm] = sbp_half(order,n,h,r)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_half(order,n,h,r)
%
% Construct SBP staggered grid operators with minimal boundary closure size 
%
% Input:
% order : Order of accuracy
% n     : Number of grid points n+1 (nodal grid) n+2 (cell-centered grid)
% h     : Grid spacing
% r     : Reflection coefficient

%
% Output
% xp,xm        : Grid vectors xp (nodal grid) xm (cell-centered grid)
% Qp,Qm, Pp,Pm : Staggered grid operators
% alpha,bet    : SAT penalty parameters
    switch order
      case 2
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_2nd(n,h);
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_4th(n,h);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_6th(n,h,r);
      case 8
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_8th(n,h,r);
      otherwise
       error('SBP staggered grid minimal operator not implemented');   
      end

xp = xp';
xm = xm';


end
