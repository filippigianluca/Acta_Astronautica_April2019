function F = mask_CUBESAT_5subsystems(x, par)

d = x(1:par.dim_d);
u = x(par.dim_d+1:end);

F = CUBESAT_5subsystems(d, u, par);

end