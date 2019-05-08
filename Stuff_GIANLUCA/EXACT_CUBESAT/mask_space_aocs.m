function P = mask_space_aocs(x, par)


d = x(1:par.dim_d);
u = x(par.dim_d+1:end);

[M,P,info] = space_aocs(d,u);
return