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

d2r = pi/180;
r2d = 180/pi;
AU2km = 149597870.66;                         % [km]
km2AU = 1/AU2km;

% Set constants
T_sid = 86164.1004;                           % [s]
mu_S = 1.327178*1e11*km2AU^3*T_sid^2;         % [AU^3/day^2]


%Celestial bodies parameters (from Novak,Vasile-2011 and from JPL ephemeris
%Horizon system - http://ssd.jpl.nasa.gov/horizons.cgi)
%Time used for the computation (also for the true anomaly) -> 01/01/2000,
%equivalent to 0 MJD2000

% Angles->radians, ang rates->rad/day, dist->AU

%Earth
body(1).name = 'Earth';
body(1).a = 1.00000011;
body(1).e = 0.01671022;
body(1).i = 0.000054*d2r;
body(1).Omega = -0.0892318; %-11.26064*d2r;
body(1).w = 1.8857; %102.94719*d2r;
body(1).M0 = 358.189*d2r;
body(1).theta0 = 358.126*d2r;
body(1).n = 9.85*1e-1*d2r;          
body(1).t0 = 0; 

%Mars
body(2).name = 'Mars';
body(2).a = 1.524;
body(2).e = 0.093;
body(2).i = 1.85*d2r;
body(2).Omega = 49.557*d2r;
body(2).w = 286.502*d2r;
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
body(3).w = 183.267*d2r;
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
body(4).w = 178.926*d2r;
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
body(5).w = 265.647*d2r;
body(5).M0 = 266.599*d2r;
body(5).theta0 = 265.325*d2r;
body(5).n = 5.969*1e-3*d2r;     
body(5).t0 = 0; 