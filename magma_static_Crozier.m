function [rho, K, c, a, b, p0,mu, CondLen]=magma_static_Crozier(z,P)
%computes magmastatic profile from code of Josh Crozier
%Crozier and Karlstrom 2022

%http://calcul-isto.cnrs-orleans.fr/thermodynamics/applications/#/script/js_compute~H2O-CO2%20systems~H2O%20and%20CO2%20solubility.js
CO2solppm_table = readmatrix('new_Values_CO2solubility_ppm_comb.csv');
H2Osolwtp_table = readmatrix('new_Values_H2Osolubility_wtp_comb.csv');


[input, out] = magma_prop_from_volatiles_variableA_noextra(P.conddip,...
    P.H2O_condtop_wtp,P.H2O_condbot_wtp,P.CO2_condtop_ppm,P.CO2_condbot_ppm,...
    P.H2O_laketop_wtp,P.H2O_lakebot_wtp,P.CO2_laketop_ppm,P.CO2_lakebot_ppm,...
    P.Hlake,P.Hcolumn,P.Rres,P.Rcondtop,P.Rcondbot,P.Rlake,...
    P.TK_condbot,P.TK_condtop,P.TK_lakebot,P.TK_laketop,...
    'true', CO2solppm_table, H2Osolwtp_table, P.dz);



rho = interp1(out.zvec_condlake,out.density_vec_CL,z);
K = interp1(out.zvec_condlake,out.bulkmod_vec_CL,z);
c = interp1(out.zvec_condlake,out.soundspeed_vec_CL,z);
mu = interp1(out.zvec_condlake,out.mu_vec_CL,z);
CO2 = interp1(out.zvec_condlake,out.CO2_vec_CL,z);
H2O = interp1(out.zvec_condlake,out.H2O_vec_CL,z);
p0 = interp1(out.zvec_condlake,out.Pvec_CL,z);
ntot = CO2+H2O; %total mass fraction gas n

bbar = interp1(out.zvec_condlake,out.dgasmassfrac_dP_vec,z);
aa=interp1(out.zvec_condlake,out.ddensity_dn_vec,z);

%compute the solubility functions a and b
b = -bbar;%gradient(ntot)./gradient(p0); % = -dn_eq/dp
a = -1./rho .*aa;%(gradient(rho)./gradient(ntot)); % = -rho^-1 * drho/dn
a(~isfinite(a))=0;b(~isfinite(b))=0;
a = 0.04*a; %reduce...

%keyboard
%smooth the results
span = 15;
rho = smooth(rho,span);
K = smooth(K,span);
c = smooth(c,span);
mu = smooth(mu,span);
p0 = smooth(p0,span);
a = smooth(a,span);
b = smooth(b,span);

CondLen = input.Hcond;

%flip the resulting vectors
rho = flipud(rho);
K = flipud(K);
c = flipud(c);
mu = flipud(mu);
p0 = flipud(p0);
a = flipud(a);
b = flipud(b);
