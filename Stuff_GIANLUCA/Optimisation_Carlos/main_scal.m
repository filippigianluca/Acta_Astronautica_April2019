% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_scal'));
% for benchmarking purposes
% nfevals_io = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];
maxfevals = [6e4 1e6 4e5 4e5 6e5 5e4];
run0 = 16;
runfi = 4;

%%
for runid= run0 + (1:runfi)
    for tc = 1:6
        % % tic
        % global nfevalglobal;
        % nfevalglobal = 0;
        
        % global history;
        % history=[];
        if tc<7
            figure()
            clf
            load(strcat('TC_',num2str(tc),'_global_matlab.mat'));
            plot(archive(:,end-1),archive(:,end),'y*');
            grid()
            drawnow
            hold on
        end
        %% initialise problem TC1
        [ problem_minmax ] = init_tc(tc);
        problem_minmax.maxnfeval = maxfevals(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% modify nfeval (for benchmarking purposes)
            % algo_outer.par.nFeValMax = nfevals_io(tc,1);
            % algo_inner.par.nFeValMax = nfevals_io(tc,2);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        savefolder = strcat('RESULTS/scal_new/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_MO_',num2str(tc),'_',num2str(runid),'.mat'));
    end
end