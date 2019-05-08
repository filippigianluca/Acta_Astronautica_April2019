% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_mso'));
% for benchmarking purposes
% nfevals_io = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];
run0 = 6;
runfi = 2;

%%
for runid= run0 + (1:runfi)
    for tc = 20
        for w = 0.0:0.1:1.0
            
            % tic
            global nfevalglobal;
            nfevalglobal = 0;
            
            global history;
            history=[];

            %% initialise problem TC1
            [ problem_minmax ] = init_tc_mso(tc);
            problem_minmax.par_objfun{1}.lambda = [w  1.0-w];

            %% initialise algorithm minmax (META-ALGO)
            [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

            %% modify nfeval (for benchmarking purposes)
                % algo_outer.par.nFeValMax = nfevals_io(tc,1);
                % algo_inner.par.nFeValMax = nfevals_io(tc,2);

            %% optimise
            [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

            %create results directory
            savefolder = strcat('RESULTS/scal/');
            mkdir(savefolder);
            save(strcat(savefolder,'testcase_results_TC1_scal_',num2str(w),'_',num2str(runid),'.mat'));
        end

    end
end