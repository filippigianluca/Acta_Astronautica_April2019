% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_so_nested_local'));
savefolder = strcat('RESULTS/SO_NESTED_LOCAL_convex/');

for runid= 2:2     % 20 runs
    for tc =13   % 13 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;
        
        global history;
        history=[];

        %% initialise problem
        [ problem_minmax ] = init_tc_so_convex(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        dmin
        fminmax
        nfevalglobal
        
        %create results directory
        mkdir(savefolder);
         save(strcat(savefolder,'testcase_results_TC_SO_convex_nested_local_',num2str(tc),'_',num2str(runid)));

    end
end