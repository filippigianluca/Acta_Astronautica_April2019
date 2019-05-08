function [N_foto2obdh, N_foto_tot, Tac] = payload(d, u, par)


G   = 6.67e-11; % constant     (NM2/kg2)
M_T = 5.97e24;  % mas earth    (kg) 
R_e = 6.37e6;   % earth radius (km)
h      % altitude, epistemic
f_foto % foto frquency, fixed
TIME 
X 


% period of the circolar orbit
T_orbit = 2*pi*((R_e+h)^3/(G*M_T))^0.5;

% access time: time when the CubeSat is visible from the ground station
Tac = T_orbit/180*acos(X);


N_loop = floor(TIME/T_orbit);
more_time = TIME/T_orbit - N_loop;



N_foto2obdh = f_foto*N_loop + f_foto*more_time*(more_time>T_orbit-Tac);
N_foto_tot  = f_foto*TIME;
end