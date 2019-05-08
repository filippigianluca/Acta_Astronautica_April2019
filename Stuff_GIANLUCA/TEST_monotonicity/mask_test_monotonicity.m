function F = mask_test_monotonicity(x,par)


d = par{1, 1}.d_fix;
u = par{1, 1}.u_fix;

u(par{1, 1}.u_fixed_Partial) = par{1, 1}.u_coupled;
u(par{1, 1}.u_opti) = x;


F = -IAC_objfun_zeno_SA(d, u, []);
return