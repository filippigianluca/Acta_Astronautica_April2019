addpath(genpath('spaceart_toolbox'))

% Equinoctial elements
x(1) = 7000e3;
x(2) = 0;
x(3) = 0;
x(4) = 0;
x(5) = 0;
x(6) = 0; % Mean longitude

date.year = 2016;
date.month = 3;
date.day = 1;
date.hour = 0;
date.minutes = 0;
date.seconds = 0;

date.MJD2000 = date2mjd2000([date.year date.month date.day ...
    date.hour date.minutes date.seconds]);

d = date.MJD2000 - date2mjd2000([date.year 1 1 ...
    date.hour date.minutes date.seconds]);




Fx = 100;
Fb = 120;
Kp = 6;
alpha_Sun = 30*pi/180;

rho = TD88(x, Fx, Fb, Kp, d, alpha_Sun);