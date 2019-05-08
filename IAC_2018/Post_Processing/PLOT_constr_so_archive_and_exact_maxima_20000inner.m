clear all; %close all; clc


%% %%%%%%%%%%%%%%%
% set parameters %
%%%%%%%%%%%%%%%%%%

plot_exact_maxima_result_minmax = 0;

nu = 400;
crosscheck_flag = 1;
num_feval_innerouter = [1000 5000 10000 20000];


% first point to be plotted from the vector results
start_point_plot = 2;

%% %%%%%%%%%%%%%%%%%%%
% end set parameters %
%%%%%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% start add folders to path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')
%warning('on','all')


% folders in UTOPIAR_RELIABLES repository
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end
ALLdir    = CURRENTPATH(1:idcs(end-3)-1);
MY_GITHUB = CURRENTPATH(1:idcs(end-2)-1);
if isunix
    Stuff_Danda = strcat(MY_GITHUB,'/utopiae-reliables/Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'/utopiae-reliables/IAC_2018');
else
    Stuff_Danda = strcat(MY_GITHUB,'\utopiae-reliables\Stuff_DANDA');
    IAC2018     = strcat(MY_GITHUB,'\utopiae-reliables\IAC_2018');
end

rmpath(genpath(ALLdir));
addpath(genpath(IAC2018));
addpath(genpath(Stuff_Danda));

% results .txt in H drive
addpath(genpath('H:\My Documents\RESULTS_UTOPIAE_RELIABLES\ACTA_ASTRONAUTICA_IAC2018'));

% tool for plot in smart-o2c repository
addpath(genpath('C:\Users\yhb17181\Documents\MY_GITHUB\smart-o2c\Other_Tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end add folders to path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%



lineStyles=linspecer(4,1);




figure
hold on

for j = num_feval_innerouter
    ind = find(num_feval_innerouter==j);
    
    
    load(strcat('nu',num2str(nu),'_IAC_2018_so_constr_nfeval_minmax_2000000nfeval_inner_outer_',num2str(j)),'minmax');

    if crosscheck_flag == 0
        Arch = minmax.output.ARCHIVE{1, 1};
    elseif crosscheck_flag == 1
        Arch = minmax.output.ARCHIVE_crosscheck{1, 1};
    end

    

    s(ind) = scatter(cell2mat(Arch(1+ start_point_plot:end,7)), cell2mat(Arch(1+start_point_plot:end,4)), 'filled','DisplayName',strcat('nfeval_{inout} = ',num2str(j)));
    s(ind).MarkerFaceColor = lineStyles(ind,:);
    
    p(ind) = plot(cell2mat(Arch(1+ start_point_plot:end,7)), cell2mat(Arch(1+ start_point_plot:end,4)),'Color',lineStyles(ind,:),'linewidth',1.5);



end

hold off

xlabel('nfeval')
ylabel('constrained Mass')
lgd = legend(s(1,:));
title(strcat('constrained Mass minmax convergence'))




if plot_exact_maxima_result_minmax
    
    for i=1:8
        
        load(strcat('exact_maximum_400nu20000inner_30000nfeval_20agents_D',num2str(i),'interior-point_nu_390'));        
        scatter(cell2mat(Arch(1+i,7)), -fval, 'filled');
        
    end
    
end