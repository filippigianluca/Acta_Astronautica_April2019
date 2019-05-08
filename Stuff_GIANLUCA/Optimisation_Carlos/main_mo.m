% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func(strcat('init_algo_mo_esteco'));
savefolder = strcat('RESULTS/MO_RR/');

%%
for runid = 1:1
    for tc = 6

        disp(strcat('MO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;
        
        
        figure(1)
        clf
        load(strcat('TC_',num2str(tc),'_global_matlab.mat'));
        plot(archive(:,end-1),archive(:,end),'y*');
        grid()
        drawnow
        hold on
        figure(1)

        %% initialise problem 
        [ problem_minmax ] = init_tc(tc);
        
        % ONLY FOR BENCHMARK REPEATABILITY PURPOSES! {{{
             init = str2func(strcat('init_algo_mo_ausa_tc',num2str(tc)));
        % }}}

        %% initialise algorithm minmax (META-ALGO).
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);
        
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);


        % figure(tc+1)
        % hold on
        % load(strcat('TC_',num2str(tc),'_global_matlab.mat'));
        % plot(archive(:,end-1),archive(:,end),'k*');
        % plot(fminmax(:,1),fminmax(:,2),'r.')

        %create results directory
        mkdir(savefolder);
%         save(strcat(savefolder,'testcase_results_TC_MO_',num2str(tc),'_',num2str(runid)));
    end
end
