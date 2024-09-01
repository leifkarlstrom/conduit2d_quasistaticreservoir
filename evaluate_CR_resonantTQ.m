%evaluate natural frequency and quality factor of conduit-reservoir mode
%based on the magmastatic model of Crozier and Karlstrom 2022

clear

%need to import background vectors from magmastatic code to do this
%justice!

%to test, assume 400 m long conduit with dz = 10 
%z is positive up (last entry is top of lake)

paramsweep = false; % set this to false for single evaluation of model, set to true for a sweep (and then define which variable)

N = 40; %number of grid points
Pm.dz = 10; %uniform spatial step used to define magmastatic profile in m
Pm.H = N*Pm.dz; %length of the conduit
Pm.rho = linspace(1500,1500,N+1);%1000*ones(1,N+1); %density as fctn of depth z kg/m3
Pm.drhodz = (Pm.rho(end)-Pm.rho(1))/Pm.H *ones(1,N+1); %density gradient 
Pm.mu = 1e1*ones(1,N+1); %viscosity as fctn of depth z in Pas
Pm.betam = 1e-9; %compressibility of magma at the base of the conduit in 1/Pa

%parameter definitions for conduit (need to be consistent with simulation!)
Pm.g = 9.8; %gravity
Pm.R0 = 5; %reference conduit radius in m 

Pm.rad = Pm.R0*ones(1,N+1); %conduit radius as fct of z

%parameter definitions for chamber
Pm.G = 5e9; %host rock shear modulus in Pa
Pm.nu = 0.25; %Poisson's ratio
Pm.betar = 3/(4*Pm.G); %chamber elastic compressiblity for sphere
Pm.Dc = 500; %diameter of spherical chamber in m

Pm.Cr = pi*Pm.Dc^3/(8*Pm.G)*(1+4/3*Pm.betam*Pm.G);%total storativity of reservoir

if paramsweep == false
    %in this case, simply evaluate the model once

%evaluate the constants needed to define properties of damped harmonic
%oscillator (the reduced order model of Crozier and Karlstrom 2022)
[c1,c2,c3,Omega] = evaluateparams(Pm);

%note! angular frequency is related to frequency by factor of 2pi!
Tvis = 2*pi/Omega;

DeltaRho = Pm.rho(end)-Pm.rho(1);

%inviscid resonant period
Tinvis = 2*pi*sqrt(Pm.H*mean(Pm.rho)/ ...
    ((Pm.rho(end)-DeltaRho)*Pm.g+pi*Pm.rad(1)^2/Pm.Cr));

%fully developed resonant period
Tfd = 2*pi./sqrt(((Pm.rho(end)-DeltaRho)*Pm.g+pi*Pm.rad(1)^2/Pm.Cr)./(Pm.H*mean(Pm.rho)) - ...
    16*mean(Pm.mu)^2/(Pm.rad(1)^4*mean(Pm.rho)^2));

lambda = 4*mean(Pm.mu)/(Pm.rad(1)^2*mean(Pm.rho));

%quality factor
Q = Omega*c1/c2;

%fully developed Quality factor
Qfd = Pm.rad(1)^2 * mean(Pm.rho)/(8*mean(Pm.mu)) * ...
    sqrt((Pm.g*(Pm.rho(end)-DeltaRho)+ pi*Pm.rad(1)^2/Pm.Cr)/(Pm.H*mean(Pm.rho)) - ...
    16*mean(Pm.mu)^2/(Pm.rad(1)^4*mean(Pm.rho)^2));

disp(['resonant period is ' num2str(Tvis) ' sec'])
disp(['Quality factor is ' num2str(Q)])

disp(['For reference, inviscid resonant period is ' num2str(Tinvis) ' sec'])
disp(['For reference, fully developed resonant period is ' num2str(Tfd) ' sec'])
disp(['For reference, fully developed Q is ' num2str(Qfd) ])

elseif paramsweep == true

%do a sweep through parameters and plot up the result 
%explore viscosity as example

MuVec = logspace(0,3,15);

for ii = 1:length(MuVec)
    disp(ii)
    Pm.mu = MuVec(ii)*ones(1,N+1); %viscosity as fctn of depth z in Pas

    [c1,c2,c3,Omega] = evaluateparams(Pm);

    %note! angular frequency is related to frequency by factor of 2pi!
    Tvis(ii) = 2*pi/Omega;

    %quality factor
    Q(ii) = Omega*c1/c2;
end

yyaxis left
plot(MuVec,Tvis,'bo-')
ylabel('Period (sec)')
yyaxis right
plot(MuVec,Q,'rs--')
ylabel('Quality factor')
xlabel('Viscosity (Pas)')
set(gca,'Xscale','log')
end

%% here are the functions we're evaluating the above script

function [c1,c2,c3,Omega] = evaluateparams(Pm)

%inertial term
c1 = Pm.rad(1).^2*trapz(Pm.dz,Pm.rho./Pm.rad.^2);
%restoring force term
c3 = -Pm.rad(1).^2*(Pm.g*(trapz(Pm.dz,Pm.drhodz./Pm.rad.^2)-(Pm.rho(end)./Pm.rad(end).^2)) - pi/Pm.Cr);

%solve for the natural angular frequency of the system
z0=0.01; %starting guess for minimization

%solve implicitly for conduit-reservoir natural angular frequency (inverse of
%period)
c2 = dampingterm(z0,Pm);

%solve implicitly for resonance including frequency dependent damping
if c2==0||isnan(c2)
    Omega = sqrt(c3/c1);
else
    Omega = fzero(@(omega) myfun(omega,Pm,c1,c3),z0);
end


c2 = dampingterm(Omega,Pm);

end


function c2 = dampingterm(omega,Pm)
%calculate viscous damping term
alpha = sqrt(omega.*Pm.rho./Pm.mu)*1i^(3/2);
integrand = (Pm.mu.*alpha)./Pm.rad.^3 .*besselj(1,Pm.rad.*alpha)./besselj(2,Pm.rad.*alpha);
c2 = 2*Pm.rad(1).^2 * real(trapz(Pm.dz,integrand));
end


function f = myfun(omega,Pm,c1,c3)

c2 = dampingterm(omega,Pm);

%throw an error if overdamped (too high viscosity)
    if c3/c1<(c2/(2*c1))^2%abs(imag(f))>0
        disp(c3/c1)
        disp((c2/(2*c1))^2)
        error('System is overdamped, no resonance')    
    end

%now set up the function to be minimized (move terms in eqn 28 all to one
%side)
f = omega - sqrt(c3/c1 - (c2/(2*c1)).^2);

end







