% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_so_new'));
savefolder = strcat('RESULTS/SO_new_2/');

for runid= 1:20     % 20 runs
    for tc = 1:13   % 13 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;
        
        global history;
        history=[];

        %% initialise problem
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        mkdir(savefolder);
         save(strcat(savefolder,'testcase_results_TC_SO_new_',num2str(tc),'_',num2str(runid)));

    end
end