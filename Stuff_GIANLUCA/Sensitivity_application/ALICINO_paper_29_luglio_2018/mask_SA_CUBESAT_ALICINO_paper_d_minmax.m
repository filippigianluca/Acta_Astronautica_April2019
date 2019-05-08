function F = mask_SA_CUBESAT_ALICINO_paper_d_minmax(x, par)


d = par.d;
u = x;

F = CUBESAT_ALICINO_paper(d,u);
return