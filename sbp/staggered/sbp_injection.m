function [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection(order,n,h,truncate)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection(order,n,h,truncate)
    switch order
      case 2
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_2nd(n,h,1);
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_4th(n,h);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_6th(n,h);
      case 8
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_8th(n,h);
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

      xp = xp';
      xm = xm';


end
