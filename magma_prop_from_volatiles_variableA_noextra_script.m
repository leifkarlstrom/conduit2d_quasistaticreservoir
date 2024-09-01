%script to calculate magmastatic properties

conddip = 90;
Rres = 750; 
% Zrescentroid = 1100-1940;
% Zlaketop = 900;
Zlaketop = 250;
% Zlakebot = 700;
Zlakebot = 0;
% Zcondbot = Zrescentroid + Rres;
% Hcolumn = Zlaketop - Zcondbot;
Hcolumn = 750; % Lake + Conduit 
% Hlake = Zlaketop - Zlakebot;
Hlake = 250; 
Rcondtop = 15;   %15
Rcondbot = 15;   %15
Rlake = 115; % For 2018: Fig 3 (Crozier&Karlstrom 2022)
Zrescentroid = Hcolumn + Rres;


% Fig 3 gives a range of ~ 1000 - 1200 deg C (Crozier&Karlstrom 2022)
Magma_temp = 1155; % one consistent temp for the system

TK_condbot = Magma_temp + 273.15; 
TK_condtop = Magma_temp + 273.15;
TK_lakebot = Magma_temp + 273.15; 
TK_laketop = Magma_temp + 273.15; 

H2O_laketop_wtp = 0.3;  %volatile contents - 1.5
H2O_lakebot_wtp = 0.3;  %volatile contents - 1.5
H2O_condtop_wtp = 0.3;  %volatile contents - 1.5 
H2O_condbot_wtp = 0.1;  %volatile contents %vary 0.1 - 0.5

wtp_CO2_to_ppm = 1/(7.0778e-05);
CO2_laketop_ppm = 1/2*H2O_laketop_wtp*wtp_CO2_to_ppm;  %ratio can change 
CO2_lakebot_ppm = 1/2*H2O_lakebot_wtp*wtp_CO2_to_ppm;
CO2_condtop_ppm = 1/2*H2O_condtop_wtp*wtp_CO2_to_ppm;
CO2_condbot_ppm = 1/2*H2O_condbot_wtp*wtp_CO2_to_ppm;

dz = 0.2; %depth increment for integration, recommend 0.4 or less

plotprofiles = false;

%http://calcul-isto.cnrs-orleans.fr/thermodynamics/applications/#/script/js_compute~H2O-CO2%20systems~H2O%20and%20CO2%20solubility.js
CO2solppm_table = readmatrix('new_Values_CO2solubility_ppm_comb.csv');
H2Osolwtp_table = readmatrix('new_Values_H2Osolubility_wtp_comb.csv');


[input_params, magmastatic_struct] = magma_prop_from_volatiles_variableA_noextra(conddip,...
    H2O_condtop_wtp,H2O_condbot_wtp,CO2_condtop_ppm,CO2_condbot_ppm,...
    H2O_laketop_wtp,H2O_lakebot_wtp,CO2_laketop_ppm,CO2_lakebot_ppm,...
    Hlake,Hcolumn,Rres,Rcondtop,Rcondbot,Rlake,...
    TK_condbot,TK_condtop,TK_lakebot,TK_laketop,...
    plotprofiles, CO2solppm_table, H2Osolwtp_table, dz);

% save("test_Temp_array.mat", "magmastatic_struct", "input_params")


