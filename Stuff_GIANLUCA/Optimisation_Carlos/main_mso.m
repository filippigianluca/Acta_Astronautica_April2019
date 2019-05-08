% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_mso'));
% for benchmarking purposes
% nfevals_io = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];
run0 = 20;
runfi = 20;

%%
for runid= run0 + (1:runfi)
    for tc = 1:7
        for obj = 1: (2 + (tc==7))

            % tic
            global nfevalglobal;
            nfevalglobal = 0;
            
            global history;
            history=[];

            tc_mso = 10*tc + obj;
            %% initialise problem TC1
            [ problem_minmax ] = init_tc_mso(tc_mso);

            %% initialise algorithm minmax (META-ALGO)
            [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

            %% modify nfeval (for benchmarking purposes)
                % algo_outer.par.nFeValMax = nfevals_io(tc,1);
                % algo_inner.par.nFeValMax = nfevals_io(tc,2);

            %% optimise
            [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

            %create results directory
            savefolder = strcat('RESULTS/MSO_new/');
            mkdir(savefolder);
            save(strcat(savefolder,'testcase_results_TC_mso_',num2str(tc_mso),'_',num2str(runid)));
        end
    end
end