% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func(strcat('init_minmarek'));
savefolder = strcat('RESULTS/EMO/');

%%
for runid=1
    for tc = 1
        disp(strcat('EMO',num2str(tc),'_',num2str(runid)))       
        
%         clf
         %figure(runid)
%         load(strcat('/home/carlos/phd/MACSMINMAX_testbench/MO/testcase_macsminmax_after_refactoring/reference_PyGMO_runs/global_fronts/matlab/TC_',num2str(tc),'_global_matlab.mat'));
%         plot(archive(:,end-1),archive(:,end),'k*');
%         grid()
%         hold on
        figure(1)
        clf
        load(strcat('TC_',num2str(tc),'_global_matlab.mat'));
        plot(archive(:,end-1),archive(:,end),'y*');
        grid()
        drawnow
        hold on
        figure(1)

        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem 
        [ problem_minmax ] = init_tc(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_EMO_',num2str(tc),'_',num2str(runid)));
    end
end

% error_runs