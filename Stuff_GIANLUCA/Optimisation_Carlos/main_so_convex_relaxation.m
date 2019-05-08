% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_convex_relaxation'));
savefolder = strcat('RESULTS/SO_convex_relaxation/');

maxnfevaltot = [5e3 1e4 1.5e4 7.5e3 1e4 1.5e4 4e4 4e2 2e3 1e3 6e3 3e5 2e6 5e4 2.5e6];
nfeval_outer = [1e2 1e2 1e2 1e2 2e2 4e2 5e2 1e2 6e1 1e2 2e2 6e2 8e3 1.5e3 4.8e3];
for runid= 1     % 20 runs
    for tc = 1% 15 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc_so_convex(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        algo_minmax.par_minmax.maxnfeval = maxnfevaltot(tc);
        algo_outer.par_au.nFeValMax = nfeval_outer(tc);
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        dmin
        fminmax
        nfevalglobal
        output.iter

        % figure(tc)
        % plot(output.convergence(:,1),output.convergence(:,2));
        % drawnow
        % hold on

        %create results directory
        mkdir(savefolder);
         save(strcat(savefolder,'testcase_results_TC_SO_convex_relaxation_',num2str(tc),'_',num2str(runid)));

    end
end