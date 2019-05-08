function F = IAC_function_tune_fix_u(d, par) 


u = par.u;
F = CUBESAT_5subsystems_MASS(d, u, par); 

% global f
% f = [f F];
return