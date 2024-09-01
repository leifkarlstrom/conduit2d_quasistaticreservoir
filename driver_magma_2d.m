function out=driver_magma_2d(nr,nz,order,InputParam,plot_solution,tot_time,bgstate)
%solver 2d magma conduit problem with viscosity and quasi-static bottom crack
% nr: number of grid points in r direction
% nz:  number of grid points in z direction is nz+1
% order: order of accuracy
% mu: viscosity
% dt  : time step
% plot_solution: true (plot), false (don't plot)
%
% material properties (stored in data structure M)
% (no specific units assumed, so either nondimensionalize
% or use SI or other self-consistent set of units)
%
% M.rho = density
% M.c = sound speed
% M.K = bulk modulus
% M.g = gravitational acceleration
% M.a = gas exsolution, -(drho/dn)/rho
% M.b = gas exsolution, -dneq/dp
% M.tau = gas exsolution time scale
% M.mu = coefficient of friction at walls of conduit

% addpath('../../matlab/sbp');
% addpath('../../matlab/sbp/staggered');
% addpath('../../matlab/io');
% addpath('../../matlab/time_integrators');
% addpath(genpath(pwd));%add all the subfolders

% by default, do not plot solution
if nargin==2, order = 2;mu = 0; plot_solution = false; end
if nargin==3, plot_solution = false; mu = 0; end
if nargin==4, plot_solution = false; end

% Domain size
%R = 30; % conduit radius
L = InputParam.Hcolumn;%700; % conduit length 
M.L = L;

R = InputParam.Rcondtop; %conduit radius

% Grid points and grid spacing
dz = L/nz; M.dz = dz;
dr = InputParam.Rcondtop/nr;%R/nr;
 

z = dz*[0:nz]'; % Construct z grid

%bottom boundary condition
M.bottom_bc = 'crack_quasi_static' ;%'p=0';%
% 'p=0': zero pressure 
% 'v=0': zero velocity
% 'crack_quasi_static': quasi-static crack/reservoir


%background state:
M.background_state = bgstate;%'CrozierKarlstrom';%'parameterized';'magma_static';
% 'constant': constant
% 'magma_static': magma static
% 'parameterized' : prescribe variation in properties

M.background_gas = 'Halemaumau_mean';%'water_shallow';%'Halemaumau_max';
%'Halemaumau_max';
%'Halemaumau_mean';
%'water_shallow';

%calculate contributions to energy
M.energy_calc=false;

% define background state parameters;
M.g = InputParam.g; 

switch M.background_state
    case 'parameterized'
        % specify bgstate coefficients

        M.rho = InputParam.rho;
        M.c = InputParam.c;

        M.K = M.rho.*M.c.^2;
        
        %parameters related to kinetic degassing model - set to 0 for now
        M.a = InputParam.a;%(1-ex);  -0.3210
        M.b = InputParam.b;%(1-ex); -5.1852e-07

        M.tau = InputParam.tau;%tau; 
        M.mu = InputParam.mu; %viscosity, constant for now

        %keyboard
    case 'magma_static'
        % magma static condition
        [M.rho, M.K, M.c, M.a, M.b, M.p0, M.ex]=magma_st(z,M.background_gas);
        M.tau = 0.05;%tau;
        M.R = 10;
        keyboard
    case 'CrozierKarlstrom'
        [M.rho, M.K, M.c, M.a, M.b, M.p0, mu, ConLen]=magma_static_Crozier(z,InputParam);
        M.mu = mean(mu(z<100)); %mean(mu);
        
        out.mu = mu;
        out.InputParam = InputParam;

        M.CondLen = ConLen; %this is Hcond from the magmastatic script

        M.tau = 0.05; %volatile diffusion timescale

end

% %build in lava lake parametrized by tanh function 
%M.R=R*ones(nz+1,1); % radius of the conduit;
M.lakedepth = M.L;%250; %depth of lake
M.lakeradius = InputParam.Rlake;%100 - R;
M.conduitradius = InputParam.Rcondtop; %assume constant

M.lakebottomflatness = .1; %scale for how sharp to make conduit radius transition

M.R=M.conduitradius*ones(nz+1,1) + (M.lakeradius-M.conduitradius)*(1+tanh(M.lakebottomflatness*(z-M.lakedepth)))/2;

M.S=pi*M.R.^2;%pi*R^2*ones(nz+1,1);% cross-section area of the conduit

%properties of wall rock
M.K_w = InputParam.K_w;% bulk modulus of wall rock
M.nu_w = InputParam.nu_w;% Poisson's ratio of wall rock
M.lambda_w = 3*M.K_w*M.nu_w/(1+M.nu_w); % lame constant of wall rock
M.G_w = 3*M.K_w*(1-2*M.nu_w)/(2+2*M.nu_w);% shear modulus of rock wall rock  

% properties of bottom crack    
if strcmp(M.bottom_bc,'crack_quasi_static')

    M.A_c = M.S(1);% cross section area at conduit bottom (do not confused with the crack surface area)
%     M.R_c = 550; % radius of the penny shape crack
%     M.W_c = 1.3; % maximum width of the penny shape crack
%     M.V_c = 4/3*pi*M.W_c*M.R_c^2; % crack volume 
%     M.K_c = 3*pi*M.K_w*(1-2*M.nu_w)/(4*M.R_c/(M.W_c/2)*(1-M.nu_w^2)); % crack stiffness. (1/V_c)*dV_c/d_p
%     M.alpha = M.A_c/M.V_c*1/((1/M.K(1)+1/M.K_c));% coupling parameters, dp_c/dt=alpha*v(0).

%sphere case
    M.R_c = InputParam.Rres;%750; % semi major radius of the spheroid
    M.W_c = InputParam.Rres;%750; %semiminor radius
    M.V_c = 4/3*pi*M.W_c*M.R_c^2; % sphere volume 
    M.K_c = 4*M.G_w/3;%1./M.V_c * M.R_c^2*pi*M.W_c/M.G_w; % sphere elastic stiffness. (1/V_c)*dV_c/d_p
    M.alpha = M.A_c/M.V_c*1/((1/M.K(1)+1/M.K_c));% coupling parameters, dp_c/dt=alpha*v(0).
%keyboard
end

% forcing at top
M.pT.A = InputParam.TopForceAmp; %30e3; % pressure perturbation amplitude (in Pa)
M.pT.T = InputParam.TopForceDur;%1; % pressure perturbation duration
M.pT.t = InputParam.TopForceCenter;%7; % pressure perturbation center time

% internal forcing (input of energy within conduit)
M.ic.A = 0e3; % Amplitude
M.ic.w = 0.2*L; % Gaussian width for initial conditions
M.ic.zc = (1/2)*L; % Gaussian center

%run the model
[A,G,F,lambda,z,r,ind,VP,W2] = magma_2d(nr,M.conduitradius,nz,L,M,order);

%keyboard
% F = F*0; % Turn off forcing

% time steps

%total simulation time
tmax=tot_time;%1*L/min(M.c);

if nargin<8 % timestep not specified
    CFL=0.5; 
    overdamp = 1; %LK added, decrease timestep below CFL for stability if needed
    dt = overdamp*CFL*dz/max(M.c);
end

% tmax = 120; %4*L/min(M.c); % Final time
nt = ceil(tmax/dt); % Number of time steps
dt = tmax/nt; % adjust dt to finish exactly at t=tmax


f = @(q,t) A*q + F*G(t);

v0 = zeros(nr*(nz+1),1);
p0 = zeros((nz+1),1);
% v0 = 0.01*ones(nr*(nz+1),1); % Constant initial velocity
% p0 = M.ic.A*exp((-0.5*(z-M.ic.zc).^2)/(M.ic.w^2));
n0 = zeros((nz+1),1);
h0 = zeros(1,1);
M_c0 = zeros(3,1) ; %moment from conduit
q0 = [v0; p0; n0; h0];

if strcmp(M.bottom_bc, 'crack_quasi_static')
    %add additional unknow p_c: crack pressure
    p_c0 = zeros(1,1);
    q0     = [v0; p0; n0; h0; p_c0];
end

q = q0;

% Storage arrays
out.v = zeros(nr*(nz+1),nt+1);
out.p = zeros((nz+1),nt+1);
out.n = zeros((nz+1),nt+1);
out.h = zeros(1       ,nt+1);
out.totalEnergy = zeros(1       ,nt+1); %total energy in the conduit
out.BGRdiss = zeros(1       ,nt+1); 
out.Viscdiss = zeros(1       ,nt+1); 
out.M_c = zeros(3, nt+1);% moment from the conduit

out.vbar = zeros(nz+1,nt+1); %cross sectional average

% save initial conditions
out.v(:,1) = v0;
out.p(:,1) = p0;
out.n(:,1) = n0;
out.h(:,1) = h0;
out.M_c(:,1) = M_c0;



out.z = z;
out.t = linspace(0,tmax,nt+1);
out.dt = dt;
out.nt = nt;
out.M = M;

if strcmp(M.bottom_bc, 'crack_quasi_static')
    out.p_c = zeros(1,1);
    out.p_c (:,1) = p_c0;
    out.M_f = zeros(3, nt+1);%moment from fracture
end

if strcmp(plot_solution,'true')
    figure('units','normalized','outerposition',[0 0 1 1]);
end

if M.energy_calc
    truncate=1;
    [rp,rm,Pp,Pm,Qp,Qm] = sbp_injection(8,nr,dr,truncate);
    Rp = spdiags(rp',0,nr+1,nr+1);
    Rm = spdiags(rm',0,nr,nr);
    [Pz_inv,Dz] = sbp_sparse(8,nz,dz);
    Pz = inv(Pz_inv);
    invPp = inv(Pp);
    em = ones([nr,1]);
end

for m=0:nt-1
        
    q = lsrk4(f,q,m*dt,dt);
    
    v = q(ind.iv);
    p = q(ind.ip);
    n = q(ind.in);
    h = q(ind.ih);
    
    % Save fields to output
    out.v(:,m+2) = v;
    out.p(:,m+2) = p;
    out.n(:,m+2) = n;
    out.h(:,m+2) = h;
    out.M_c(:, m+2) = [M.lambda_w+2*M.G_w,M.lambda_w+2*M.G_w,M.lambda_w]'/M.G_w*...
                                2*pi*dz*sum([0.5;ones(nz-1,1);0.5].*p.*M.R.^2);

    %width averaged velocity
    vb=W2*v;
    %Vbar(:,m+2) = vb;
    out.vbar(:,m+2) = vb;

    if strcmp(M.bottom_bc,'crack_quasi_static')
        p_c = q(ind.ip_c);
        out.p_c(:, m+2) = p_c;
        out.M_f(:,m+2) = 16*(1-M.nu_w^2)*M.R_c^3/(9*M.K_w*(1-2*M.nu_w))...
          *p_c*[M.lambda_w, M.lambda_w, M.lambda_w+2*M.G_w]'; 
    end
    
    if M.energy_calc %calculate the total energy in the conduit at each timestep

% %width averaged velocity
% vb=W2*v;
% %Vbar(:,m+2) = vb;
% out.vbar(:,m+2) = vb;

[out.BGRdiss(m+2),out.Viscdiss(m+2),KE,PE,NE,GE] = dissipation(v,p,n,h,M,R,Pz,mu,Rp,Pp,invPp,Qp,VP,Rm,Pm,nz,em);  

out.totalEnergy(m+2) = KE+PE+NE+GE;

if mu > 0
    out.viscoustime = R.^2.*M.rho./mu;
end

% if m+2>2 && TotE(m+2)>TotE(m+1) 
%     maxE = TotE(m+2); %find max energy
% end
% if m+2>2 && TotE(m+2) < maxE/exp(1)
%     out.endInd = m+2;
%     keyboard
% end
%plot(Viscdiss(:,m+2));drawnow
%keyboard
    end

    if strcmp(plot_solution,'true')

        Rlake = dr*[0:(nr)*round(InputParam.Rlake/InputParam.Rcondtop)-1];    

        if mod(m,50)==0
            vmat = reshape(v,[nr,nz+1]);

            condarray = NaN(nr*round(InputParam.Rlake/InputParam.Rcondtop),nz+1);


            %condarray(1:nr,1:ConLen/dz)=vmat(:,1:ConLen/dz);

            for II = 1:length(vmat(1,:))%1:(length(vmat(1,:))-(ConLen/dz+1))
                Rfac = round(out.M.R(II)/InputParam.Rcondtop);
                
                %dB = repmat(vmat(:,ConLen/dz+II)',(InputParam.Rlake/InputParam.Rcondtop),1);
                dB = repmat(vmat(:,II)',Rfac,1);

                %conrow = NaN(size(condarray(:,1)));
%                 if max(abs(dB))>0
%                 keyboard
%                 end
                condarray(1:nr*Rfac,II)=reshape(dB,[1,nr*Rfac]);

                %condarray(:,ConLen/dz+II) = reshape(dB,[1,nr*(InputParam.Rlake/InputParam.Rcondtop)]);
                %condarray(:,II) = reshape(dB,[1,nr*Rfac]);
            end
            



            subplot(1,3,1)
            pcolor(Rlake,z,condarray'); shading flat
            %pcolor(r,z,vmat'), shading flat
%             colorbar
            title('Velocity')
            xlabel('r (m)'), ylabel('z (m)')
            clim([-.3 .3])
            
            subplot(1,3,2)
            plot(p,z), xlim([-7e4,7e4])
            ylim([0, L])
            title('Pressure')
            xlabel('p'), ylabel('z (m)')
            
            subplot(1,3,3)
            plot(n,z)
            xlim([-5e-3, 5e-3])
            ylim([0, L])
            title('Exsolved Gas Mass Fraction')
            xlabel('n'), ylabel('z (m)')
            %pause(1e-9)
            drawnow();
           
        end
    end
    if mod(m,1000)==0
        disp(['time = ' num2str(m*dt)])
    end
end


%keyboard
if strcmp(plot_solution,'true')
    figure;
    subplot(4,1,1)
    plot(out.t,out.h);xlabel('time');ylabel('h');
    subplot(4,1,2)
    plot(out.t,out.p(25,:));xlabel('time');ylabel(['Conduit pressure at z=' num2str(out.z(25)) ' m'])
    subplot(4,1,3)
    plot(out.t,out.p_c);xlabel('time');ylabel('P_c')

    %now plot spectrum
L = length(out.p(2,:));                                           % length of signal
Fs = 1/out.dt;                                              % Make Up Sampling Frequency & Units (Hz)
Fn = Fs/2;                                              % Nyquist Frequency
FTs = fft(out.p(10,:) - mean(out.p(10,:)))/L;              % Subtract Mean, look right above chamber
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
Iv = 1:numel(Fv);                                       % Index Vector

%hold off
semilogx(Fv, abs(FTs(Iv))*2)
grid
xlabel('Frequency (Hz)')
ylabel('Amplitude in conduit')


    figure;
    pcolor(out.t,out.z,out.p);shading flat
    xlabel('time (s)'); ylabel('depth')
    cmap
    caxis([-5e4,5e4]); colorbar


end

end
