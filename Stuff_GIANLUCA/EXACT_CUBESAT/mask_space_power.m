function M = mask_space_power(x, par)


d = x(1:par.dim_d);
u = x(par.dim_d + 1:end);

[M,P,info] = space_power(d,u);
return