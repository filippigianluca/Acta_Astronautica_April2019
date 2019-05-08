function F = mask_SA_CUBESAT_ALICINO_paper(x, par)


d = x(1:par.dim_d);
u = x(par.dim_d+1:end);

F = CUBESAT_ALICINO_paper(d,u);
return