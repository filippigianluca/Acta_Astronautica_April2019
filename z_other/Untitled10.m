load('fix_d_u_par');
for i=1:1000
    
    par.fix.time = i;
    [F] = CUBESAT_5subsystems(d, u, par);
    out(i) = F;
end