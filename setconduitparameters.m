function [MSIn] = setconduitparameters(bgstate,nz)

% % load parameter file
% MSIn = load("NonInverted_Params_04062018_T202927.mat");
            % fixed (see paper) 290

MSIn.g = 9.8; %gravity m/s

MSIn.Hcolumn = 500;  % lake + conduit in m
MSIn.Rcondtop = 10;  % conduit radius in m
MSIn.Rlake = 10;    % lake radius in m

dz = MSIn.Hcolumn/nz; M.dz = dz;
z = dz*[0:nz]'; % Construct z grid 

%elastic properties of wall rock    
MSIn.K_w =10e9;% bulk modulus of wall rock (Pa)
MSIn.nu_w = 0.25;% Poisson's ratio of wall rock

MSIn.Rres = 200; %reservoir radius in m

MSIn.TopForceAmp = 30e3; %pressure perturbation amplitude (in Pa)
MSIn.TopForceDur = 1; % pressure perturbation duration
MSIn.TopForceCenter = 5; % pressure perturbation center time



if strcmp(bgstate,'parameterized')   

    option = 'step';

    %different functional forms for the transition in properties
    MSIn.c0 = 400; %reference sound speed (m/s)
    MSIn.rho0 = 1500; %reference density (kg/m3)

    switch option
        case 'const'
            MSIn.c = MSIn.c0 * ones(length(z),1);
            MSIn.rho = MSIn.rho0 * ones(length(z),1);
        case 'step'
            Lex = 0.5*MSIn.Hcolumn; % transition depth
            %ex = heaviside(1-z/Lex); % abrupt transition
            ex = 1-tanh((z/Lex).^10); % smoothed transition, varies between 0 and 1
            MSIn.rho = 500+1000*ex; 
            MSIn.c = MSIn.c0+400*ex; %sound speed   
        case 'exp'
        % parameterization in Liang et al 2019a
        % in this case you specify delta rho
        % note: Liang neglects inertial contribution of Lake in CR formula,
        % which is not justified for straight conduit (need large lake)
   
            MSIn.rhoL = 500;
            MSIn.alpha = MSIn.Hcolumn/(log(MSIn.rho0)-log(MSIn.rhoL));
            MSIn.rho = MSIn.rhoL*exp(( MSIn.Hcolumn-z)/MSIn.alpha); %slightly different than Liang since z increases down for us
            MSIn.c = MSIn.c0 * ones(length(z),1); %Liang chooses constant sound speed
    end

        MSIn.K = MSIn.rho.*MSIn.c.^2;
        
        %parameters related to kinetic degassing model in KD2016 - set to 0 for now
        MSIn.a = 0*ones(length(z),1);%(1-ex);  -0.3210
        MSIn.b = 0*ones(length(z),1);%(1-ex); -5.1852e-07

        MSIn.tau = .05;%tau volatile diffusion time (from KD16) in sec; 
        MSIn.mu = 0; %viscosity, constant for now

end



