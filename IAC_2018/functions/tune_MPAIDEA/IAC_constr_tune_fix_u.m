function [F, Feq] = IAC_constr_tune_fix_u(d, par) 


u = par.u;
F = CUBESAT_5subsystems_RES(d, u, par) + par.fix.nu;

Feq = [];
return