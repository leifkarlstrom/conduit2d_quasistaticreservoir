function [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_8nd(n,h,x)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_injection_8nd(n,h,x)
% n : number of grid points  (n+1) nodal grid (n+2) cell-centered grid
% h : grid spacing

% Unknown coefficients
x = [...
   1.086496652252120;
  -0.057499088529918;
   0.005314209550716;
   0.012075871596527;
   0.107508816058426;
   1.141192616998850;
  -0.015219347808173;
   0.017375374161306;
  -0.013153917729437;
   1.014984881252215;
   ];
pm7 = x(10);
qm74 = x(8);
qm01 = x(1);
qm07 = x(3);
qm05 = x(2);
qm72 = x(7);
qm62 = x(4);
qm64 = x(5);
qm75 = x(9);
qm67 = x(6);


% Coefficients determined such that the SBP property is satisfied
qm44 = -2675*pm7/64 - 25*qm01/8 - 75*qm05/8 + 175*qm07/4 + 10*qm64 + 40*qm74 + 534420319/12386304;
qp44 = -1415*pm7/32 - 25*qm01/4 - 75*qm05/4 + 175*qm07/2 + 10*qm64 + 45*qm74 + 85868021/1769472;
qm31 = -101*pm7 - 10*qm01 + 6*qm62 - 2*qm64 - 14*qm67 + 30*qm72 - 18*qm74 - 15*qm75 + 8076005953/61931520;
qp71 = -qm07;
pm4 = -35*pm7 + 17740819/483840;
qm56 = 7215*pm7/64 - qm01/24 + 5*qm05/24 + 35*qm07/12 + qm62/3 - 2*qm64/3 + 56*qm67/3 - 5*qm72/3 + 6*qm74 + 40*qm75/3 - 375186975923/2786918400;
qp54 = 2359*pm7/8 + 10*qm05 + 2*qm62 - 6*qm64 + 42*qm67 + 45*qm75 - 35694471229/103219200;
qm76 = -367*pm7/120 + qm72/9 - 2*qm74/5 - 8*qm75/9 + 1009567/322560;
qm06 = 17*pm7/576 + qm01/24 - 5*qm05/24 - 35*qm07/12 - 9138299/139345920;
pp3 = -35*pm7 + 830952317/23224320;
qm43 = -2185*pm7/24 + 25*qm01/3 + 25*qm05/3 - 175*qm07/3 - 20*qm62/3 - 20*qm64/3 - 70*qm67/3 - 200*qm72/9 - 40*qm74 - 200*qm75/9 + 61843386637/557383680;
qp27 = -qm62;
qm66 = -6769*pm7/360 - qm62/15 + 2*qm64/15 - 56*qm67/15 + 31764529/1433600;
qm63 = -1127*pm7/72 - 2*qm62/3 - 2*qm64/3 - 7*qm67/3 + 5306543/286720;
qp60 = 0;
qp00 = -1;
qp61 = -17*pm7/576 - qm01/24 + 5*qm05/24 + 35*qm07/12 + 9138299/139345920;
qm53 = 435*pm7/8 - 5*qm01/3 - 5*qm05/3 + 35*qm07/3 + 10*qm62/3 + 10*qm64/3 + 35*qm67/3 + 25*qm72/3 + 15*qm74 + 25*qm75/3 - 931992349/13934592;
qm04 = -31*pm7/64 - 5*qm01/8 - 15*qm05/8 + 35*qm07/4 + 16663957/15482880;
qm73 = 3*pm7/2 - 5*qm72/9 - qm74 - 5*qm75/9 - 48533/32256;
qm23 = -271*pm7/4 + 50*qm01/3 + 50*qm05/3 - 350*qm07/3 - 10*qm62/3 - 10*qm64/3 - 35*qm67/3 - 40*qm72/3 - 24*qm74 - 40*qm75/3 + 36686200759/557383680;
qm46 = -15579*pm7/64 + 5*qm01/24 - 25*qm05/24 - 175*qm07/12 - 2*qm62/3 + 4*qm64/3 - 112*qm67/3 + 40*qm72/9 - 16*qm74 - 320*qm75/9 + 806956984979/2786918400;
qm12 = 1777*pm7/64 + 75*qm01/8 + 25*qm05/8 - 105*qm07/4 - qm62 - 5*qm72 - 2480055887/66355200;
qp10 = 0;
qp33 = 271*pm7/4 - 50*qm01/3 - 50*qm05/3 + 350*qm07/3 + 10*qm62/3 + 10*qm64/3 + 35*qm67/3 + 40*qm72/3 + 24*qm74 + 40*qm75/3 - 36686200759/557383680;
qp73 = -6187*pm7/168 - 10*qm07 - 5*qm67 + 8*qm72/7 - 24*qm74/7 - 40*qm75/7 + 841232540273/19508428800;
qp38 = -3*pm7/2 + 5*qm72/9 + qm74 + 5*qm75/9 + 48533/32256;
qp55 = -2177*pm7/8 - 5*qm05 - 2*qm62 + 6*qm64 - 42*qm67 - 40*qm75 + 299458197427/928972800;
qp53 = -1183*pm7/8 - 10*qm05 - qm62 + 3*qm64 - 21*qm67 - 24*qm75 + 32207984797/185794560;
qp48 = -qm74;
qp52 = 1183*pm7/40 + 5*qm05 + qm62/5 - 3*qm64/5 + 21*qm67/5 + 5*qm75 - 32207336797/928972800;
qp24 = -293*pm7/32 - 75*qm01/4 - 25*qm05/4 + 105*qm07/2 + 10*qm62 + 45*qm72 + 9380843393/309657600;
qm70 = 89*pm7/168 + 10*qm72/63 - qm74/7 - 8*qm75/63 - 1198679/2257920;
qp32 = -8*pm7/3 + 25*qm01/3 + 25*qm05/3 - 175*qm07/3 - 2*qm62/3 - 2*qm64/3 - 7*qm67/3 - 25*qm72/9 - 5*qm74 - 25*qm75/9 - 1647131291/557383680;
qp25 = 265*pm7/64 + 75*qm01/8 + 25*qm05/8 - 105*qm07/4 - 10*qm62 - 40*qm72 - 2736330311/185794560;
qm17 = -107*pm7/14 - 5*qm07 - qm67 + 5*qm72/21 - 5*qm74/7 - 25*qm75/21 + 57906662941/6502809600;
pm0 = -pm7 + 537547/241920;
pp2 = 21*pm7 - 779227069/38707200;
qm37 = -1947*pm7/28 - 10*qm07 - 10*qm67 + 15*qm72/7 - 45*qm74/7 - 75*qm75/7 + 177905115601/2167603200;
qp70 = 0;
pp4 = 35*pm7 - 155921497/4644864;
qp17 = -371*pm7/40 + 3*qm62/5 - qm64/5 - 7*qm67/5 + 15753779/1433600;
qp64 = -74803*pm7/288 + 5*qm01/12 - 25*qm05/12 - 175*qm07/6 - 2*qm62/3 + 4*qm64/3 - 112*qm67/3 + 5*qm72 - 18*qm74 - 40*qm75 + 170386623871/557383680;
pm2 = -21*pm7 + 4484093/193536;
qp16 = 101*pm7/3 + qm01 - 3*qm62 + qm64 + 7*qm67 - 10*qm72 + 6*qm74 + 5*qm75 - 60482587699/1393459200;
qp65 = 15579*pm7/64 - 5*qm01/24 + 25*qm05/24 + 175*qm07/12 + 2*qm62/3 - 4*qm64/3 + 112*qm67/3 - 40*qm72/9 + 16*qm74 + 320*qm75/9 - 806956984979/2786918400;
qp41 = 31*pm7/64 + 5*qm01/8 + 15*qm05/8 - 35*qm07/4 - 16663957/15482880;
qm77 = 129*pm7/56 - qm72/21 + qm74/7 + 5*qm75/21 - 526499/150528;
qp51 = -qm05;
qm42 = -265*pm7/64 - 75*qm01/8 - 25*qm05/8 + 105*qm07/4 + 10*qm62 + 40*qm72 + 2736330311/185794560;
pp1 = -7*pm7 + 966361469/116121600;
qp57 = -287*pm7/10 - qm62/5 + 3*qm64/5 - 21*qm67/5 + 36317567/1075200;
qp67 = 6769*pm7/360 + qm62/15 - 2*qm64/15 + 56*qm67/15 - 31764529/1433600;
qp20 = 0;
qp36 = -435*pm7/8 + 5*qm01/3 + 5*qm05/3 - 35*qm07/3 - 10*qm62/3 - 10*qm64/3 - 35*qm67/3 - 25*qm72/3 - 15*qm74 - 25*qm75/3 + 931992349/13934592;
qm40 = -38461*pm7/1344 - 25*qm01/24 + 5*qm05/24 - 25*qm07/12 + 4*qm62/3 - 2*qm64/3 - 16*qm67/3 + 400*qm72/63 - 40*qm74/7 - 320*qm75/63 + 711150296867/19508428800;
qp21 = 109*pm7/64 + 15*qm01/8 + 5*qm05/8 - 21*qm07/4 - 58592623/15482880;
qm65 = 287*pm7/10 + qm62/5 - 3*qm64/5 + 21*qm67/5 - 36317567/1075200;
qp05 = 38461*pm7/1344 + 25*qm01/24 - 5*qm05/24 + 25*qm07/12 - 4*qm62/3 + 2*qm64/3 + 16*qm67/3 - 400*qm72/63 + 40*qm74/7 + 320*qm75/63 - 711150296867/19508428800;
qm33 = 8245*pm7/72 - 50*qm01/3 - 50*qm05/3 + 350*qm07/3 + 20*qm62/3 + 20*qm64/3 + 70*qm67/3 + 25*qm72 + 45*qm74 + 25*qm75 - 10093449097/79626240;
qm61 = 371*pm7/40 - 3*qm62/5 + qm64/5 + 7*qm67/5 - 15753779/1433600;
qp62 = -26429*pm7/960 + 5*qm01/24 - 25*qm05/24 - 175*qm07/12 - qm62/15 + 2*qm64/15 - 56*qm67/15 + 5*qm72/9 - 2*qm74 - 40*qm75/9 + 44557081159/1393459200;
qp42 = -155*pm7/64 - 25*qm01/8 - 75*qm05/8 + 175*qm07/4 + qm64 + 5*qm74 + 16663957/3096576;
qm26 = -21559*pm7/160 + 5*qm01/12 - 25*qm05/12 - 175*qm07/6 - qm62/3 + 2*qm64/3 - 56*qm67/3 + 8*qm72/3 - 48*qm74/5 - 64*qm75/3 + 439011863021/2786918400;
qp58 = -qm75;
qp43 = 605*pm7/96 + 25*qm01/4 + 75*qm05/4 - 175*qm07/2 - 5*qm64 - 24*qm74 - 1355726401/111476736;
qm71 = -51*pm7/40 - 2*qm72/3 + 2*qm74/5 + qm75/3 + 137401/107520;
qp13 = -2477*pm7/30 - 10*qm01 + 3*qm62 - qm64 - 7*qm67 + 16*qm72 - 48*qm74/5 - 8*qm75 + 288949094201/2786918400;
qp23 = 1021*pm7/32 + 75*qm01/4 + 25*qm05/4 - 105*qm07/2 - 5*qm62 - 24*qm72 - 48402433973/928972800;
qm11 = -1201*pm7/40 - 5*qm01 + 3*qm62/5 - qm64/5 - 7*qm67/5 + 10*qm72/3 - 2*qm74 - 5*qm75/3 + 34229931131/928972800;
qp76 = 4153*pm7/168 + qm07 + 5*qm67 - 5*qm72/7 + 15*qm74/7 + 25*qm75/7 - 150135268049/4877107200;
qm50 = 5715*pm7/448 + 5*qm01/24 - qm05/24 + 5*qm07/12 - 2*qm62/3 + qm64/3 + 8*qm67/3 - 50*qm72/21 + 15*qm74/7 + 40*qm75/21 - 317410231637/19508428800;
qm60 = -1267*pm7/360 + 2*qm62/15 - qm64/15 - 8*qm67/15 + 5981571/1433600;
pp0 = pm7 - 78498587/116121600;
qm57 = -4153*pm7/168 - qm07 - 5*qm67 + 5*qm72/7 - 15*qm74/7 - 25*qm75/7 + 150135268049/4877107200;
qp31 = -25*pm7/18 - 5*qm01/3 - 5*qm05/3 + 35*qm07/3 + 2687735/870912;
qp26 = -53*pm7/64 - 15*qm01/8 - 5*qm05/8 + 21*qm07/4 + 5*qm62 + 15*qm72 + 2736330311/928972800;
qp30 = 0;
qp06 = -5715*pm7/448 - 5*qm01/24 + qm05/24 - 5*qm07/12 + 2*qm62/3 - qm64/3 - 8*qm67/3 + 50*qm72/21 - 15*qm74/7 - 40*qm75/21 + 317410231637/19508428800;
qm34 = 1415*pm7/32 + 25*qm01/4 + 75*qm05/4 - 175*qm07/2 - 10*qm64 - 45*qm74 - 85868021/1769472;
qp07 = 1267*pm7/360 - 2*qm62/15 + qm64/15 + 8*qm67/15 - 5981571/1433600;
qm35 = -2359*pm7/8 - 10*qm05 - 2*qm62 + 6*qm64 - 42*qm67 - 45*qm75 + 35694471229/103219200;
qm25 = 1183*pm7/8 + 10*qm05 + qm62 - 3*qm64 + 21*qm67 + 24*qm75 - 32207984797/185794560;
pp6 = 7*pm7 - 681977789/116121600;
qp56 = 497*pm7/4 + qm05 + qm62 - 3*qm64 + 21*qm67 + 15*qm75 - 69209136721/464486400;
qm24 = -605*pm7/96 - 25*qm01/4 - 75*qm05/4 + 175*qm07/2 + 5*qm64 + 24*qm74 + 1355726401/111476736;
qp11 = -qm01;
qp66 = -7215*pm7/64 + qm01/24 - 5*qm05/24 - 35*qm07/12 - qm62/3 + 2*qm64/3 - 56*qm67/3 + 5*qm72/3 - 6*qm74 - 40*qm75/3 + 375186975923/2786918400;
qm13 = 8*pm7/3 - 25*qm01/3 - 25*qm05/3 + 175*qm07/3 + 2*qm62/3 + 2*qm64/3 + 7*qm67/3 + 25*qm72/9 + 5*qm74 + 25*qm75/9 + 1647131291/557383680;
qp46 = -373*pm7/192 - 5*qm01/8 - 15*qm05/8 + 35*qm07/4 + 5*qm64 + 15*qm74 + 1379509937/557383680;
qp15 = -593*pm7/8 - 5*qm01 + 6*qm62 - 2*qm64 - 14*qm67 + 80*qm72/3 - 16*qm74 - 40*qm75/3 + 12907725239/132710400;
qp47 = -qm64;
qp72 = 107*pm7/14 + 5*qm07 + qm67 - 5*qm72/21 + 5*qm74/7 + 25*qm75/21 - 57906662941/6502809600;
qm14 = 155*pm7/64 + 25*qm01/8 + 75*qm05/8 - 175*qm07/4 - qm64 - 5*qm74 - 16663957/3096576;
qm54 = 373*pm7/192 + 5*qm01/8 + 15*qm05/8 - 35*qm07/4 - 5*qm64 - 15*qm74 - 1379509937/557383680;
pp7 = -pm7 + 229947227/116121600;
qp68 = 367*pm7/120 - qm72/9 + 2*qm74/5 + 8*qm75/9 - 1009567/322560;
qm00 = 443*pm7/576 - 5*qm01/24 + qm05/24 - 5*qm07/12 - 238133321/139345920;
qp08 = -89*pm7/168 - 10*qm72/63 + qm74/7 + 8*qm75/63 + 1198679/2257920;
qm36 = 74803*pm7/288 - 5*qm01/12 + 25*qm05/12 + 175*qm07/6 + 2*qm62/3 - 4*qm64/3 + 112*qm67/3 - 5*qm72 + 18*qm74 + 40*qm75 - 170386623871/557383680;
qp12 = 1201*pm7/40 + 5*qm01 - 3*qm62/5 + qm64/5 + 7*qm67/5 - 10*qm72/3 + 2*qm74 + 5*qm75/3 - 34229931131/928972800;
qm45 = 2177*pm7/8 + 5*qm05 + 2*qm62 - 6*qm64 + 42*qm67 + 40*qm75 - 299458197427/928972800;
qp35 = 2185*pm7/24 - 25*qm01/3 - 25*qm05/3 + 175*qm07/3 + 20*qm62/3 + 20*qm64/3 + 70*qm67/3 + 200*qm72/9 + 40*qm74 + 200*qm75/9 - 61843386637/557383680;
qp45 = 2675*pm7/64 + 25*qm01/8 + 75*qm05/8 - 175*qm07/4 - 10*qm64 - 40*qm74 - 534420319/12386304;
qp22 = -1777*pm7/64 - 75*qm01/8 - 25*qm05/8 + 105*qm07/4 + qm62 + 5*qm72 + 2480055887/66355200;
qp02 = -15363*pm7/2240 - 25*qm01/24 + 5*qm05/24 - 25*qm07/12 + 2*qm62/15 - qm64/15 - 8*qm67/15 + 50*qm72/63 - 5*qm74/7 - 40*qm75/63 + 21079214299/2438553600;
qm30 = 76183*pm7/2016 + 25*qm01/12 - 5*qm05/12 + 25*qm07/6 - 4*qm62/3 + 2*qm64/3 + 16*qm67/3 - 50*qm72/7 + 45*qm74/7 + 40*qm75/7 - 917228048351/19508428800;
qp18 = 51*pm7/40 + 2*qm72/3 - 2*qm74/5 - qm75/3 - 137401/107520;
pm5 = 21*pm7 - 19565701/967680;
pp5 = -21*pm7 + 838678909/38707200;
qm51 = -101*pm7/3 - qm01 + 3*qm62 - qm64 - 7*qm67 + 10*qm72 - 6*qm74 - 5*qm75 + 60482587699/1393459200;
qp03 = 5951*pm7/224 + 25*qm01/12 - 5*qm05/12 + 25*qm07/6 - 2*qm62/3 + qm64/3 + 8*qm67/3 - 80*qm72/21 + 24*qm74/7 + 64*qm75/21 - 126982260409/3901685760;
qm21 = 2477*pm7/30 + 10*qm01 - 3*qm62 + qm64 + 7*qm67 - 16*qm72 + 48*qm74/5 + 8*qm75 - 288949094201/2786918400;
qp40 = 0;
qm10 = 15363*pm7/2240 + 25*qm01/24 - 5*qm05/24 + 25*qm07/12 - 2*qm62/15 + qm64/15 + 8*qm67/15 - 50*qm72/63 + 5*qm74/7 + 40*qm75/63 - 21079214299/2438553600;
pm3 = 35*pm7 - 2428439/69120;
qm20 = -5951*pm7/224 - 25*qm01/12 + 5*qm05/12 - 25*qm07/6 + 2*qm62/3 - qm64/3 - 8*qm67/3 + 80*qm72/21 - 24*qm74/7 - 64*qm75/21 + 126982260409/3901685760;
qp75 = -3515*pm7/56 - 5*qm07 - 10*qm67 + 40*qm72/21 - 40*qm74/7 - 200*qm75/21 + 97821536333/1300561920;
qp74 = 1947*pm7/28 + 10*qm07 + 10*qm67 - 15*qm72/7 + 45*qm74/7 + 75*qm75/7 - 177905115601/2167603200;
qm03 = 25*pm7/18 + 5*qm01/3 + 5*qm05/3 - 35*qm07/3 - 2687735/870912;
qp01 = -443*pm7/576 + 5*qm01/24 - qm05/24 + 5*qm07/12 + 238133321/139345920;
qm52 = 53*pm7/64 + 15*qm01/8 + 5*qm05/8 - 21*qm07/4 - 5*qm62 - 15*qm72 - 2736330311/928972800;
qp50 = 0;
qp78 = -129*pm7/56 + qm72/21 - qm74/7 - 5*qm75/21 + 526499/150528;
qm02 = -109*pm7/64 - 15*qm01/8 - 5*qm05/8 + 21*qm07/4 + 58592623/15482880;
qm47 = 3515*pm7/56 + 5*qm07 + 10*qm67 - 40*qm72/21 + 40*qm74/7 + 200*qm75/21 - 97821536333/1300561920;
qm22 = -1021*pm7/32 - 75*qm01/4 - 25*qm05/4 + 105*qm07/2 + 5*qm62 + 24*qm72 + 48402433973/928972800;
qm15 = -1183*pm7/40 - 5*qm05 - qm62/5 + 3*qm64/5 - 21*qm67/5 - 5*qm75 + 32207336797/928972800;
qp37 = 1127*pm7/72 + 2*qm62/3 + 2*qm64/3 + 7*qm67/3 - 5306543/286720;
pm6 = -7*pm7 + 287831/35840;
qp28 = -qm72;
qp04 = -76183*pm7/2016 - 25*qm01/12 + 5*qm05/12 - 25*qm07/6 + 4*qm62/3 - 2*qm64/3 - 16*qm67/3 + 50*qm72/7 - 45*qm74/7 - 40*qm75/7 + 917228048351/19508428800;
qm32 = 293*pm7/32 + 75*qm01/4 + 25*qm05/4 - 105*qm07/2 - 10*qm62 - 45*qm72 - 9380843393/309657600;
qm41 = 593*pm7/8 + 5*qm01 - 6*qm62 + 2*qm64 + 14*qm67 - 80*qm72/3 + 16*qm74 + 40*qm75/3 - 12907725239/132710400;
qm27 = 6187*pm7/168 + 10*qm07 + 5*qm67 - 8*qm72/7 + 24*qm74/7 + 40*qm75/7 - 841232540273/19508428800;
qm16 = 26429*pm7/960 - 5*qm01/24 + 25*qm05/24 + 175*qm07/12 + qm62/15 - 2*qm64/15 + 56*qm67/15 - 5*qm72/9 + 2*qm74 + 40*qm75/9 - 44557081159/1393459200;
qp14 = 101*pm7 + 10*qm01 - 6*qm62 + 2*qm64 + 14*qm67 - 30*qm72 + 18*qm74 + 15*qm75 - 8076005953/61931520;
pm1 = 7*pm7 - 6518441/967680;
qp63 = 21559*pm7/160 - 5*qm01/12 + 25*qm05/12 + 175*qm07/6 + qm62/3 - 2*qm64/3 + 56*qm67/3 - 8*qm72/3 + 48*qm74/5 + 64*qm75/3 - 439011863021/2786918400;
qm55 = -497*pm7/4 - qm05 - qm62 + 3*qm64 - 21*qm67 - 15*qm75 + 69209136721/464486400;
qp77 = -qm67;
qp34 = -8245*pm7/72 + 50*qm01/3 + 50*qm05/3 - 350*qm07/3 - 20*qm62/3 - 20*qm64/3 - 70*qm67/3 - 25*qm72 - 45*qm74 - 25*qm75 + 10093449097/79626240;



% Number of coefficients
b = 8;

% Q+ and Q-, top-left corner
QpL = [...
qp00, qp01, qp02, qp03, qp04, qp05, qp06, qp07, qp08;
 qp10, qp11, qp12, qp13, qp14, qp15, qp16, qp17, qp18;
 qp20, qp21, qp22, qp23, qp24, qp25, qp26, qp27, qp28;
 qp30, qp31, qp32, qp33, qp34, qp35, qp36, qp37, qp38;
 qp40, qp41, qp42, qp43, qp44, qp45, qp46, qp47, qp48;
 qp50, qp51, qp52, qp53, qp54, qp55, qp56, qp57, qp58;
 qp60, qp61, qp62, qp63, qp64, qp65, qp66, qp67, qp68;
 qp70, qp71, qp72, qp73, qp74, qp75, qp76, qp77, qp78
];
QmL = [...
0, 0, 0, 0, 0, 0, 0, 0;
 qm00, qm01, qm02, qm03, qm04, qm05, qm06, qm07;
 qm10, qm11, qm12, qm13, qm14, qm15, qm16, qm17;
 qm20, qm21, qm22, qm23, qm24, qm25, qm26, qm27;
 qm30, qm31, qm32, qm33, qm34, qm35, qm36, qm37;
 qm40, qm41, qm42, qm43, qm44, qm45, qm46, qm47;
 qm50, qm51, qm52, qm53, qm54, qm55, qm56, qm57;
 qm60, qm61, qm62, qm63, qm64, qm65, qm66, qm67;
 qm70, qm71, qm72, qm73, qm74, qm75, qm76, qm77
];

% Q+ and Q-
w = b; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qp(end,:) = [];
Qm(:,end) = [];

% Add SBP boundary closures
Qp(1:b,1:b+1) = QpL;
Qp(end-b+1:end,end-b:end) = -fliplr(flipud(QpL));
Qm(1:b+1,1:b) = QmL;
Qm(end-b:end,end-b+1:end) = -fliplr(flipud(QmL));

% P+ and P-
Pp = ones(n+1,1);
Pm = ones(n+2,1);

Pp(1:b) = [pp0,  pp1,  pp2,  pp3,  pp4,  pp5,  pp6,  pp7]; 
Pp(end-b+1:end) = Pp(b:-1:1);
Pm(1:b+1) = [0,  pm0,  pm1,  pm2,  pm3,  pm4,  pm5,  pm6,  pm7];
Pm(end-b:end) = Pm(b+1:-1:1);
Pp = spdiags(Pp,0,n+1,n+1);
Pm = spdiags(Pm,0,n+2,n+2);

Pp = h*Pp;
Pm = h*Pm;

% nodal and cell-centered grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  


% Test operators
test = true;
if test
for j=0:b/2
  disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
  num2str(max(abs(Qp*xm.^j-j*Pp*xp.^max([j-1,0]))))]);
  disp([ 'Dm, j = ' num2str(j) ' Error max = '...
  num2str(max(abs(Qm*xp.^j-j*Pm*xm.^max([j-1,0]))))]);
end  
end
