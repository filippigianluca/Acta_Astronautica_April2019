function [M_batteries, battery, DOD] = mask_batteries_sizing(d, e_sizing_vector, eff)


V_bus = d(1);
E_cap = d(2);


eta.system.bdr = eff;
eta.system.lcl = eff;
eta.system.harn = eff;
eta.system.batt = eff;



e_sizing_vector.in(1) = E_cap;

e_start(1) = E_cap;
for i=2:length(e_sizing_vector.in)
    e_start(i) = min(E_cap, e_start(i-1) - e_sizing_vector.out(i-1) + e_sizing_vector.in(i));
end

E_req = max(e_start - e_sizing_vector.out);
[M_batteries, battery, DOD] = batteries_sizing(E_req, E_cap, V_bus, eta);


return
