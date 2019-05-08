% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_convex_ei_outer'));
savefolder = strcat('RESULTS/SO_convex_ei_outer/');

maxnfevaltot = [2e3 2e3 3e3 2e3 2e3 4e3 2e4 1e2 5e2 2e3 3e3 1e4 2e4 2e4 2e4];
nfeval_outer = [3e3 3e3 3e3 3e3 3e3 3e3 3e3 3e3 3e3 3e3 3e3 1e4 1e4 5e3 1e4];
for runid= 1:20     % 20 runs
    for tc = 1:15% 15 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc_so_convex(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        algo_minmax.par_minmax.maxnfeval = maxnfevaltot(tc);
        algo_outer.par_sa.nFeValMax = nfeval_outer(tc);
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        dmin
        fminmax
        nfevalglobal
        output.iter

        figure(tc)
        plot(output.convergence(:,1),output.convergence(:,2));
        drawnow
        hold on

        %create results directory
        mkdir(savefolder);
         save(strcat(savefolder,'testcase_results_TC_SO_convex_ei_outer_',num2str(tc),'_',num2str(runid)));

    end
end