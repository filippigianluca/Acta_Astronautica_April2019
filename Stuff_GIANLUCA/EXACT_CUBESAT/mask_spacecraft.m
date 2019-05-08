function M = mask_spacecraft(x, par)


d_aocs = x(1:par.dim_d_aocs);
u_aocs = x(par.dim_d+1 : par.dim_d+par.dim_u_aocs);


[M_aocs, P_aocs, ~] = space_aocs(d_aocs, u_aocs);


d_ttc = x(par.dim_d_aocs +1 : par.dim_d_ttc);
u_ttc = x(par.dim_d+par.dim_u_aocs + 1 : par.dim_d + par.dim_u_ttc);


[M_ttc, P_ttc, ~] = space_ttc(d_ttc, u_ttc);


d_power = x(par.dim_d_ttc + 1 : par.dim_d);

u_power(1) = 16 + P_aocs + P_ttc;
u_power(2) = 16 + P_aocs + P_ttc;
u_power(3:12) = x(par.dim_d+par.dim_u_ttc + 1 : end);


[M_power, ~, ~] = space_power(d_power, u_power);



M = M_aocs + M_ttc + M_power;

return