clear 
close all

addpath(genpath('../'))


% Distance unit DU - Earth Radius [km]
DU = 6378.136;

% Define function to fit
% f = 'exponential_atm_model_DUTU';
f = 'MSIS_atm_model';

% Results of f are in kg/m^3
flag_log = 0;

% Degree of the polynomial
n = 4;

% Name to save results
save_filename = '../Chebyshev_110km_1000km_MSIS_test.mat';


%% 110-125

% Interval [km]
a_km =110;
b_km =125;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,1) = coefficients1;

%% 125-150

% Interval [km]
a_km =125;
b_km =150;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,2) = coefficients1;

%% 150-250

% Interval [km]
a_km =150;
b_km =250;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,3) = coefficients1;


%% 250-350

% Interval [km]
a_km =250;
b_km =350;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,4) = coefficients1;

%% 350-500

% Interval [km]
a_km =350;
b_km =500;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,5) = coefficients1;


%% 500-700

% Interval [km]
a_km =500;
b_km  = 700;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,6) = coefficients1;


%% 700-1000

% Interval [km]
a_km = 700;
b_km  = 1000;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,7) = coefficients1;


%% 1000-1500

% Interval [km]
a_km = 1000;
b_km  = 1500;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,8) = coefficients1;

%% 1500-2000

% Interval [km]
a_km = 1500;
b_km  = 2000;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,9) = coefficients1;



%% 2000-2500

% Interval [km]
a_km = 2000;
b_km  = 2500;

a = a_km/DU;
b = b_km/DU;

coefficients1 = Chebyshev_interpolation(f,n,a,b,flag_log);
coeff(:,10) = coefficients1;



%% At higher altitude: no drag!
coeff(:,11) = zeros(5,1);

% close all

save(save_filename, 'coeff')
