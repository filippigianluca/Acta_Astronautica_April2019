% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_sc'));
savefolder = strcat('RESULTS/SC_check_max/');

%% initialisations
load('d_lhs');
n_runs = 10;

[ problem_minmax ] = init_problem_sc;
[ problem_max_u  ] = build_metaproblem_macsminmax_inner(problem_minmax);
problem_max_u.par_objfun.objective = 1;

algo_max = init_algo_sc_max(problem_minmax);
 
%%
for d_i = 2
    for run_i = 1:10
        disp(strcat('SC_',num2str(d_i),'_',num2str(run_i)))
        % d = d_lhs(d_i,:)
        % d = 1e2*[0.100000000000000   0.900000000000000   0.070000000000000   0.004613197932412   0.004066392398145   0.001510415753850   1.450000000000000  0.035999596861513   0.010000000000000   0.000305482028773];
        d = 1e2*[0.100000000000000   0.900000000000000   0.070000000000000   0.004613197932412   0.004066392398145   0.003000000000000   1.450000000000000  0.035999596861513   0.010000000000000   0.000305482028773];
        
        d = (d-problem_minmax.lb_d')./(problem_minmax.ub_d'-problem_minmax.lb_d');
        problem_max_u.par_objfun.d = d;


        global nfevalglobal;
        nfevalglobal = 0;
        
        [ umax, fmax , exitflag , output] = algo_max{1}.optimise(problem_max_u,algo_max{1}.par);

        %create results directory
        mkdir(savefolder);
        save(strcat(savefolder,'SC_check_max_',num2str(d_i),'_',num2str(run_i)));
    end
end