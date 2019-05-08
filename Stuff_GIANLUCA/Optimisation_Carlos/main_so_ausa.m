% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_ausa'));
savefolder = strcat('RESULTS/SO_RBSAC/');

maxnfevaltot = [2e3 6e3 1e4 3e3 3e3 5.5e3 2.5e4 3e2 2.5e3 2e3 4e3 1.5e3 1e4]; % these are more or less the fevals with nested local psi
for runid= 1:1     % 20 runs
    for tc = 1:13   % 13 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        algo_minmax.par_minmax.maxnfeval = maxnfevaltot(tc)*100;

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        dmin
        fminmax
        nfevalglobal
        iter = output.iter

        figure(tc)
        plot(output.convergence(:,1),output.convergence(:,2));
        drawnow
        hold on

        %create results directory
        mkdir(savefolder);
         save(strcat(savefolder,'testcase_results_TC_SO_RBSAC_',num2str(tc),'_',num2str(runid)));

    end
end