%surface displacement timeseries from pipe and crack
%using output from magma2d
%clear
N=200;
%out=driver_magma_2d(35,2^8,8,0,15,1);

out.M.R=20;

model='open_pipe'; % model 'open_pipe' or 'close pipe'
N_st=1; % number of station
N_z=N+1; % number of grid points in z direction
R=out.M.R(1); % conduit radius
%c1 =0 ; % depth 1
%c2 = 20*R; % depth 2
%P=1e4; % magnitude pressure
pipe_loc=[0 0]'; % conduit location
x=0;%zeros(1,N_st); % x coordinates (East)
y = 1000;%linspace(R,10*R,N_st); %y coordinates (North)
station_loc = [x; y]; % station location
z = out.z;%linspace(c1,c2,N_z)';% depth of each grid point
mu = .4e10; % shear modulus
nu = 0.25; % Poisson's ratio

%for the crack
%x=zeros(1,1000);% East
%y=linspace(0,1e3,1000); %North
%station_loc = [x;y]; % station location
Hc=max(out.z); % depth of crack

strike=0; % strike
dip=0;  % dip

clear Upipe UcrackF t Ud
disp('calculating displacements')
ind=1;
for ii=1:20:length(out.t)
    t(ind)=out.t(ii);
    p=out.p(:,ii);%P*ones(1,N_z)'; % pressure along the depth    
    
Upipe(:,ind) = disp_pipe( pipe_loc, station_loc, flipud(p), z, R, mu, nu, model);

if strcmp(out.M.bottom_bc,'crack_quasi_static')
    pc=out.p_c(ii); % crack hydrostatic pressure
    Rc=out.M.R_c; % crack radius
Wc=out.M.W_c; % width of the crack
[U1_fialko, U2_fialko, U3_fialko] = disp_crack(Hc, Rc, Wc, strike, dip, station_loc, pc, mu, nu, 'Fialko2001'); 
%[U1_sun, U2_sun, U3_sun] = disp_crack(Hc, Rc, Wc, strike, dip, station_loc, pc, mu, nu, 'Sun1969'); 

UcrackF(1,ind)=U1_fialko;
UcrackF(2,ind)=U2_fialko;
UcrackF(3,ind)=U3_fialko;

% UcrackS(1,ind)=U1_sun;
% UcrackS(2,ind)=U2_sun;
% UcrackS(3,ind)=U3_sun;

Ud(:,ind) = Upipe(:,ind)+UcrackF(:,ind);
else
Ud(:,ind) = Upipe(:,ind);    
end
%disp(ind)
 ind=ind+1;
end

% 
figure(1)
%subplot(1,3,3)
plot(t,Ud(2:3,:))


