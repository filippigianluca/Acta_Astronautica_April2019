function F = constr_mask_CUBESAT_5subsystems_RES(d, u, par) 

F = CUBESAT_5subsystems_RES(d, u, par) + par.fix.nu;


end