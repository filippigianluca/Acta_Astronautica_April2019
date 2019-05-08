function [F, Feq] = IAC_constr_tune_fix_d(u, par) 


d = par.d;
F = CUBESAT_5subsystems_RES(d, u, par) + par.fix.nu; % CUBESAT_5subsystems_RES(d, u, par) evaluate "- RES"

Feq = [];
return