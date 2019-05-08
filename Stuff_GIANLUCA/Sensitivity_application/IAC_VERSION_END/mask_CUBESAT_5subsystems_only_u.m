function F = mask_CUBESAT_5subsystems_only_u(x, par)

d = par.d;
u = x;

F = CUBESAT_5subsystems(d, u, par);

end