% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
%init = str2func(strcat('init_algo_mo_noncon_1'));
savefolder = strcat('RESULTS/TRANSFER/');

%%
init = str2func('init_algo_transfer');
for runid = 1:20

    disp(strcat('TRANSFER_',num2str(runid)))
    global nfevalglobal;
    nfevalglobal = 0;
    
    
    figure(1)
    clf
    load('dv_transfer_new.mat','dv_min');
    plot(log10(dv_min),'y.-');
    grid()
    drawnow
    hold on
    figure(1)

    %% initialise problem 
    [ problem_minmax ] = init_problem_transfer();

    %% initialise algorithm minmax (META-ALGO).
    [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);
    
    %% optimise
    [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

    %create results directory
    mkdir(savefolder);
     save(strcat(savefolder,'TRANSFER_PIROGOV_SOUTEN_',num2str(runid)));

    figure(runid+1)
    hold on
    if runid==1
        plot(log10(dv_min),'k.-');
    end
    plot(fminmax(:,1),fminmax(:,2),'r.')

end
