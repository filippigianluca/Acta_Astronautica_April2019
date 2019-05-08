% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_so_sstl_corners'));
savefolder = strcat('RESULTS/SSTL/');
dmax = 0.5;
objective = 1;


global nfevalglobal;
nfevalglobal = 0;

global f_history;
% global theta_history;
f_history=[];
% theta_history =[];

opt_count = 0;

%% initialise problem
[ problem_0 ] = init_problem_sstl();
% problem_0.split_coordinates = [];
% problem_0.provisional_fmax = nan;

problem_list_next = [problem_0];
% n_int_orig = cellfun('size',problem_0.lb_u{objective},2);
%% build structure for approx belief curve
lb_u = problem_0.lb_u{objective}; 
ub_u = problem_0.ub_u{objective};
bpa = problem_0.bpa{objective};
dim_u = problem_0.dim_u;
n_int = cellfun('size',lb_u,2);
n_fe_tot = prod(n_int)
theta_list = [(1:n_fe_tot)',nan(n_fe_tot,dim_u+1), ones(n_fe_tot,1)];
umax_list = nan(n_fe_tot,dim_u);
for i=1:n_fe_tot
    if(mod(i,100) == 0);
        i
    end
    problem_fe = problem_0;
    pos = nfe2pos(i,n_int);
    theta_list(i,2:dim_u+1) = pos;
    for j = 1:length(pos)
        theta_list(i,end) = theta_list(i,end) * bpa{j}(pos(j));
        problem_fe.lb_u{objective}{j} = problem_0.lb_u{objective}{j}(1,pos(j));
        problem_fe.ub_u{objective}{j} = problem_0.ub_u{objective}{j}(1,pos(j));
    end
    [ ~, ~, algo_inner ] = init(problem_fe);
    f_history=[];
    problem_max_u = build_metaproblem_macsminmax_inner(problem_fe);
    problem_max_u.par_objfun.objective = objective;
    problem_max_u.par_objfun.d = (dmax-problem_fe.lb_d')./(problem_fe.ub_d'-problem_fe.lb_d');
    [ umax, fmax , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
    theta_list(i,end-1) = -fmax;
    umax_list(i,:) = map_affine(umax,problem_max_u.par_objfun.map_u_info{objective});
end


save('SSTL_bel_exact_mu05');
plot_belief(theta_list(:,end-1),theta_list(:,end));
title(num2str('Exact Belief SSTL'));