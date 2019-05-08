% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func(strcat('init_algo_nested'));

%%
for runid = 24:4:40
    for tc = 3
        disp(strcat(num2str(tc),'_',num2str(runid)))

        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem 
        [ problem_minmax ] = init_tc(tc);

        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);
        algo_inner{1}.par.nFeValMax = 603;
        algo_inner{2}.par.nFeValMax = 47463;
        algo_outer.par.maxnfeval = 120;
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        figure(tc+1)
        hold on
        load(strcat('/home/carlos/phd/MACSMINMAX_testbench/MO/testcase_macsminmax_after_refactoring/reference_PyGMO_runs/global_fronts/matlab/TC_',num2str(tc),'_global_matlab.mat'));
        plot(archive(:,end-1),archive(:,end),'k*');
        plot(fminmax(:,1),fminmax(:,2),'r.')
        drawnow
        
        %create results directory
        savefolder = strcat('RESULTS/MO_nested_cec/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_MO_',num2str(tc),'_',num2str(runid)));
    end
end

% error_runs