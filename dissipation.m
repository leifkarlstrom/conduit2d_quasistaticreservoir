function [BGRdiss, DISS,KE,PE,NE,GE] = dissipation(v,p,n,h,M,R,H,mu,Rp,Pp,invPp,Qp,VP,Rm,Pm,nz,em)
%using output of magma 2D, calculate contributions to dissipation and
%energy integrals

%dz = z(end)/nz;
%dr = R/nr;

%Rad=dr*[1:nr];

%ez = ones([nz+1,1]);
%ez = spalloc(nz+1,1,1); ez(1)= 1; ez(end)= 1;
%Ir = speye(nr);
%Iz = speye(nz+1);

a_div_b=zeros(size(M.b));
a_div_b(M.b>0)=M.a(M.b>0)./(M.b(M.b>0));
%dVdr=D1r*v;

%Viscousdiss=mu*kron(ez'*H,Ir)*(D1r*v).^2; %summed over depth but not radius
%SumVD=mu*kron(H,Ir)*(D1r*v).^2;
BGRdiss = real(ctranspose(pi*R.^2*a_div_b./M.tau.*(n + M.b.*p))*H*((n + M.b.*p)));

%Viscousdiss=2*pi*sum(Viscousdiss.*Rad');

%Viscdiss(m+2) = real(ctranspose(mu.*vb)*Pz*vb);
rho = spdiags(M.rho,0,nz+1,nz+1);
%keyboard
KE = (1/2)*2*pi*v'*kron(H*rho,Rm*Pm)*v;

%(1/2)*2*pi*transpose(v.*kron(H*M.rho,Rm*Pm*em))*v;
%real(ctranspose(0.5.*v.*VP)*v);        
PE = real(ctranspose(0.5.*(R^2*pi).*(1./M.K).*p)*H*p);
NE = real(ctranspose(0.5.*(R^2*pi).*(a_div_b).*n)*H*n);
GE = real(0.5.*(R^2*pi).*M.rho(end).*M.g.*h.^2);

%truncate = 1;
% [rp,rm,Pp,Pm,Qp,Qm] = sbp_injection(8,nr,dr,truncate);
% Rp = spdiags(rp',0,nr+1,nr+1);
% Rm = spdiags(rm',0,nr,nr);
% [Pz_inv,Dz] = sbp_sparse(8,nz,dz);
% Pz = inv(Pz_inv);

% Let's build the 1D dissipation operator     int (dv/dr)^2 r dr \approx D^T*P*R*D^T
DISS1D = 2*pi*mu*(invPp*Qp)'*Pp*Rp*invPp*Qp;

% Now let's make it 2D by integrating in the vertical direction
DISS2D = kron(H,DISS1D);

% Finally compute the total dissipation
DISS = v'*DISS2D*v; 

% if max(abs(v(:)))>1e0
%     keyboard
% end


end
