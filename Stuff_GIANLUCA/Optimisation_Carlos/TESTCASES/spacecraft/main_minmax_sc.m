% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_sc'));
savefolder = strcat('RESULTS/SC/');


%%
for runid= 1:12     % 20 runs

    disp(strcat('SC_',num2str(runid)))
    global nfevalglobal;
    nfevalglobal = 0;
    
    global history;
    history=[];

    %% initialise problem
    [ problem_minmax ] = init_problem_sc;

    %% initialise algorithm minmax (META-ALGO)
    [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

    %% optimise
    [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

    %create results directory
    mkdir(savefolder);
    save(strcat(savefolder,'SC_NEW3_',num2str(runid)));

end