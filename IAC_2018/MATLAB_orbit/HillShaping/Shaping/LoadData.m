%File with constants and parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------%
% Conversion factors
% Declare as global variables
global d2r;
global r2d;
global AU2km;
global km2AU;
global T_sid;
global mu_S;
global AU2m;
global Day2s;
global mu_T; 

d2r = pi/180;
r2d = 180/pi;
AU2km = 149597870.66;                         % [km]
km2AU = 1/AU2km;
AU2m = AU2km*1000 ; 

% Set constants
T_sid = 86164.1004;    % [s]
Day2s = T_sid;          % [s]

mu_S = 1.327178*1e11*km2AU^3*T_sid^2;         % [AU^3/day^2]
% mu_S = 1.327178*1e11;                       % [km^3/s^2]

mu_T_day = 3.986004418*1e5*T_sid^2 ; % [km^3/day^2]
mu_T = 3.986004418*1e5*(1e-4)^3 ; % [C^3/s^2] C = 10000 km


%Celestial bodies parameters (from Novak,Vasile-2011 and from JPL ephemeris
%Horizon system - http://ssd.jpl.nasa.gov/horizons.cgi)
%Time used for the computation (also for the true anomaly) -> 01/01/2000,
%equivalent to 0 MJD2000

% Angles->radians, ang rates->rad/day, dist->AU

%Earth
body(1).name = 'Earth';
body(1).a = 1.00000011 ;
% body(1).a = 1.00000011*AU2km;
body(1).e = 0.01671022;
body(1).i = 0.000054*d2r;
body(1).Omega = -0.0892318; %-11.26064*d2r;
body(1).om = 1.8857; %102.94719*d2r;
body(1).M0 = 358.189*d2r;
body(1).theta0 = 358.126*d2r;
body(1).n = 9.85*1e-1*d2r;          
body(1).t0 = 0; 

%Mars
body(2).name = 'Mars';
body(2).a = 1.524 ;
% body(2).a = 1.524*AU2km;
body(2).e = 0.093;
body(2).i = 0.85*d2r;
body(2).Omega = 49.557*d2r;
body(2).om = 286.502*d2r;
body(2).M0 = 19.095*d2r;
body(2).theta0 = 23.02*d2r;
body(2).n = 5.24*1e-1*d2r;   
body(2).t0 = 0; 


%1989ML
body(3).name = '1989ML';
body(3).a = 1.272;
body(3).e = 0.137;
body(3).i = 4.378*d2r;
body(3).Omega = 104.411*d2r;
body(3).om = 183.267*d2r;
body(3).M0 = 108.286*d2r;
body(3).theta0 = 122.256*d2r;
body(3).n = 6.865*1e-1*d2r;
body(3).t0 = 0; 

%Tempel-1
body(4).name = 'Tempel-1';
body(4).a = 3.124;
body(4).e = 0.517;
body(4).i = 10.527*d2r;
body(4).Omega = 68.933*d2r;
body(4).om = 178.926*d2r;
body(4).M0 = 359.71*d2r;
body(4).theta0 = 358.93*d2r;
body(4).n = 1.789*1e-1*d2r;
body(4).t0 = 0; 

%Neptune
body(5).name = 'Neptune';
body(5).a = 30.104;
body(5).e = 0.011;
body(5).i = 1.768*d2r;
body(5).Omega = 131.794*d2r;
body(5).om = 265.647*d2r;
body(5).M0 = 266.599*d2r;
body(5).theta0 = 265.325*d2r;
body(5).n = 5.969*1e-3*d2r;     
body(5).t0 = 0; 

% Earth Orbit 1
body(6).name = 'Orbit1';
body(6).a = 1.0000; %km
body(6).e = 0;
body(6).i = 51*d2r;
body(6).Omega = 0;
body(6).om = 0;
body(6).M0 = 0;
body(6).theta0 = 0;
body(6).n = sqrt( mu_T_day / body(6).a^3 );     
body(6).t0 = 0; 

% Earth Orbit 2
body(7).name = 'Orbit2';
body(7).a = 2.4200; %km
body(7).e = 0;
body(7).i = 56*d2r;
body(7).Omega = 150*d2r;
body(7).om = 0;
body(7).M0 = 0;
body(7).theta0 = 0;
body(7).n = sqrt( mu_T_day / body(7).a^3 );     
body(7).t0 = 0; 

%67P/Churyumov-Gerasimenko
body(8).name = '67P/Churyumov-Gerasimenko';
body(8).a = 3.464805313920435;
body(8).e = .6414365761974745;
body(8).i = 7.0452981812567*d2r;
body(8).Omega = 50.0846669914027*d2r;
body(8).om = 12.8421019463821*d2r;
body(8).M0 = 52.6900038374776*d2r;
body(8).theta0 = 265.325*d2r;
body(8).n = .152821978766721*d2r;
body(8).t0 = 6048; 

%Adonis
body(9).name = 'Adonis';
body(9).a = 1.8746121;
body(9).e = 0.76418957;
body(9).i = 1.32856*d2r;
body(9).Omega = 349.79505*d2r;
body(9).om = 43.27876*d2r;
body(9).M0 = 158.3683699*d2r;
body(9).theta0 = 265.325*d2r;
body(9).n = sqrt(mu_S / 1.8746121^3); % Only if mu_S [AU^3/day^2]
body(9).t0 = 6056; 

%Oljato            
body(10).name = 'Oljato';
body(10).a = 2.1749020;
body(10).e = 0.71309565;
body(10).i = 2.52234*d2r;
body(10).Omega = 75.00214*d2r;
body(10).om = 98.24443*d2r;
body(10).M0 = 132.8209997*d2r;
body(10).theta0 = 265.325*d2r;
body(10).n = sqrt(mu_S / 2.1749020^3); % Only if mu_S [AU^3/day^2]
body(10).t0 = 6056; 

%Epona
body(11).name = 'Epona';
body(11).a = 1.5049421;
body(11).e = 0.70223722;
body(11).i = 29.20695*d2r;
body(11).Omega = 235.50764*d2r;
body(11).om = 49.68757*d2r;
body(11).M0 = 64.2048909*d2r;
body(11).theta0 = 265.325*d2r;
body(11).n = sqrt(mu_S / 1.5049421^3); % Only if mu_S [AU^3/day^2]
body(11).t0 = 6056; 

%1974 MA
body(12).name = '1974 MA';
body(12).a = 1.7854383;
body(12).e = 0.76214174;
body(12).i = 38.06832*d2r;
body(12).Omega = 302.29037*d2r;
body(12).om = 126.91312*d2r;
body(12).M0 = 198.6814522*d2r;
body(12).theta0 = 265.325*d2r;
body(12).n = sqrt(mu_S / 1.7854383^3); % Only if mu_S [AU^3/day^2]
body(12).t0 = 6056; 

%Heracles
body(13).name = 'Heracles';
body(13).a = 1.8335239;
body(13).e = 0.77237915;
body(13).i = 9.03288*d2r;
body(13).Omega = 309.54510*d2r;
body(13).om = 227.73788*d2r;
body(13).M0 = 294.9388469*d2r;
body(13).theta0 = 265.325*d2r;
body(13).n = sqrt(mu_S / 1.8335239^3); % Only if mu_S [AU^3/day^2]
body(13).t0 = 6056; 

%Hermes
body(14).name = 'Cuno';
body(14).a = 1.9824651;
body(14).e = 0.63439136;
body(14).i = 6.70539*d2r;
body(14).Omega = 294.91792*d2r;
body(14).om = 236.30701*d2r;
body(14).M0 = 197.0602191*d2r;
body(14).theta0 = 265.325*d2r;
body(14).n = sqrt(mu_S / 1.9824651^3); % Only if mu_S [AU^3/day^2]
body(14).t0 = 6056; 