function [F] = IAC_function_MASS_tune_fix_d(u, par) 

d = par.d;
F = - CUBESAT_5subsystems_MASS(d, u, par); 

return