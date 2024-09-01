%script for running 2d conduit code with quasistatic reservoir at bottom

%based on code from https://bitbucket.org/ooreilly/magmagpu/src/master/
%associated paper: https://pangea.stanford.edu/~edunham/publications/Prochnow_etal_axisymmetric-waves_CF17.pdf

clear 
close all

%add subdirectories to path
folder = fileparts(which('run_script.m'));
addpath(genpath(folder));

%solver 2d magma conduit problem with (constant) viscosity and quasi-static bottom crack
nr = 2^4; %number of grid points in r direction % originally 2^4
nz = 2^8;  %number of grid points in z direction is nz+1 %orignally 2^8
order = 6; %order of accuracy

tot_time = 200; %total simulation time in sec 

% specify background state 
bgstate = 'parameterized';

%call script to either load parameter file or specify them, based on BGstate
[MSIn] = setconduitparameters(bgstate,nz);

tic

%run time-dependent code
out=driver_magma_2d(nr,nz,order,MSIn,'false',tot_time,bgstate);%(nr,nz,order,mu,R,tau);

toc
%keyboard
%LakeDensityAvg = mean(out.M.rho(out.z>out.M.CondLen));
%disp(['Avg lake density is ' num2str(LakeDensityAvg)])

%% post processing

% Compute the Fourier Transform of 
%chamber Pc
[Fv,FTs,Iv,spectrum] = compute_fft(out.t,out.p_c,out.dt);
%in conduit 
[Fv2,FTs2,Iv2,spectrum2] = compute_fft(out.t,out.p(3,:),out.dt);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1
periods = 1./Fv;
out.periods = periods;
out.spectrum = spectrum;
out.spectrum2 = spectrum2;

plotsolutionfields(out)

% Identify peaks in spectrum
[peaks, peak_locs]=findpeaks(spectrum,'MinPeakProminence',1.5);
%[peaks, peak_locs] = findpeaks(spectrum, 'MinPeakHeight', 0.5); % Adjust threshold as needed
disp(['peaks in spectrum ' num2str(periods(peak_locs)) ' sec'])

out.Ct = out.M.V_c*((1/out.M.K(1)+1/out.M.K_c)); %sphere storativity used at bottom of conduit

[CRout] = CR_rom_crozierkarlstrom(out);

%this implements eqn S36 in Crozier/Karlstrom 
% CRmode_constR = 2*pi*sqrt(out.M.L*0.5*(out.M.rho(end)+out.M.rho(1))/...
%     (out.M.g*(out.M.rho(end)-(out.M.rho(end)-out.M.rho(1))) + pi*out.M.R(1)^2/out.Ct));

%but version with mean of full density seems more accurate
% CRmode_constR = 2*pi*sqrt(out.M.L*mean(out.M.rho)/...
%     (out.M.g*(out.M.rho(end)-(out.M.rho(end)-out.M.rho(1))) + pi*out.M.R(1)^2/out.Ct));


OrganPipe_OO = 2*out.M.L/mean(out.M.c);

%CRmode_constR
disp(['Organ Pipe open-open 1st mode based on mean c is ' num2str(OrganPipe_OO) ' sec'])
