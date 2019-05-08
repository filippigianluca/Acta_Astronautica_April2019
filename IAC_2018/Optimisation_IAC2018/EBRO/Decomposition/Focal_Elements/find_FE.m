function FE_minmax = find_FE(n_obj, u, problem)

bpa = 1;


for i = 1:problem.dim_u
    
    for j = 1:length(problem.lb_u{n_obj}{i})
        lbj = problem.lb_u{n_obj}{i}(j);
        ubj = problem.ub_u{n_obj}{i}(j);
        if u(i) >= lbj && u(i) <= ubj
            lb(i) = lbj;
            ub(i) = ubj;
            position_FE(i) = j;
            bpa = bpa*problem.bpa{n_obj}{i}(j);
            break
        end
    end    
end




FE_minmax.lb = lb;
FE_minmax.ub = ub;
FE_minmax.bpa = bpa;
FE_minmax.position_FE  = position_FE;
return