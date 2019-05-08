clear
close all



%--------------------------------------------------------------------------
CURRENTPATH=pwd;
if isunix
    idcs   = strfind(CURRENTPATH,'/');
else
    idcs   = strfind(CURRENTPATH,'\');
end
newdir = CURRENTPATH(1:idcs(end-1)-1); %go back two folders
addpath(genpath(newdir));
%--------------------------------------------------------------------------





% rmpath(genpath('C:\Users\yhb17181\Documents\LOCAL_REPOSITORY_SMARTO2C\smart-o2c\Problems\EBRO\WorkInProgress'));
% addpath('C:\Users\yhb17181\Documents\LOCAL_REPOSITORY_SMARTO2C\smart-o2c\Problems\EBRO\WorkInProgress\CUBESAT');
% addpath(genpath('C:\Users\yhb17181\Documents\LOCAL_REPOSITORY_SMARTO2C\smart-o2c\Problems\EBRO\WorkInProgress\EXACT_CUBESAT'));


flag = 'all';

if strcmp(flag, 'aocs')
    %--------------------------------------------------------------------------
    % AOCS
    objFun = @mask_space_aocs;
    [lb_d, ub_d, lb_u, ub_u] = bound_AOCS();
    %--------------------------------------------------------------------------
elseif strcmp(flag, 'ttc')
    % TTC
    objFun = @mask_space_ttc;
    [lb_d, ub_d, lb_u, ub_u] = bound_TTC();
    %--------------------------------------------------------------------------
elseif strcmp(flag, 'power')
    % POWER
    objFun = @mask_space_power;
    [lb_d, ub_d, lb_u, ub_u] = bound_POWER();
    %--------------------------------------------------------------------------
elseif strcmp(flag, 'all')
    % ALL SYSTEM
    objFun = @mask_spacecraft;
    [lb_d, ub_d, lb_u, ub_u, par] = bound_SPACECRAFT();
    %--------------------------------------------------------------------------
end


N_x= length(lb_d) + length(lb_u);


lb = [lb_d lb_u];
ub = [ub_d ub_u];
par.dim_d = length(lb_d);





%%%%%%%%%%%%%%%%%%%%%%
distributionType='Uniform';
order=1;
printFlag=1;
%%%%%%%%%%%%%%%%%%%%%%

tic

samples=[1e4];   %[1e2 1e3]

N_runs=size(samples,2);
SV=zeros(N_x,N_runs);
ST=zeros(N_x,N_runs);
SVij=zeros(N_x,N_x,N_runs);
STij=zeros(N_x,N_x,N_runs);
A=zeros(samples(end),N_x,N_runs);
yA=zeros(samples(end),N_runs);

for i=1:N_runs
    [SV(:,i), ST(:,i), SVij(:,:,i), STij(:,:,i), A(1:samples(i),:,i), yA(1:samples(i),i)]=sensitivityAnalysis(objFun,N_x,order,lb,ub,distributionType,samples(i), par);
end

pathSA=[CURRENTPATH,'/SA_results'];
pathSAsub=[pathSA,'/',flag];
if ~(7==exist(pathSA))
    mkdir(pathSA)
end
if ~(7==exist(pathSAsub))
    mkdir(pathSAsub)
end


rawDATA=[A(:,:,end),yA(:,end)];
if strcmp(flag, 'aocs')
%--------------------------------------------------------------------------
% AOCS
filename = ['AOCS_',num2str(samples(end)),'.dat'];
figfile = fullfile(pathSAsub, filename);
csvwrite(figfile,rawDATA)
%--------------------------------------------------------------------------
elseif strcmp(flag, 'ttc')
% TTC
filename = ['TTC_',num2str(samples(end)),'.dat'];
figfile = fullfile(pathSAsub, filename);
csvwrite(figfile,rawDATA)
%--------------------------------------------------------------------------
elseif strcmp(flag, 'power')
% POWER
filename = ['POWER_',num2str(samples(end)),'.dat'];
figfile = fullfile(pathSAsub, filename);
csvwrite(figfile,rawDATA)
%--------------------------------------------------------------------------
elseif strcmp(flag, 'all')
% SpaceCraft
filename = ['ALL_',num2str(samples(end)),'.dat'];
figfile = fullfile(pathSAsub, filename);
csvwrite(figfile,rawDATA)
%--------------------------------------------------------------------------
end





figure
scatter(1:size(SV,1), SV(:,end),'filled')
title('SV')
figfile = fullfile(pathSAsub, 'SV.png');
saveas(gcf,figfile)

figure
scatter(1:size(SV,1), ST(:,end),'filled')
title('ST')
figfile = fullfile(pathSAsub, 'ST.png');
saveas(gcf,figfile)


[rank,index]=sort(ST(:,N_runs),'descend');

%identify relevent variables
relevantVars=1;
while ~(sum(rank(1:relevantVars))/sum(rank)>0.95)
relevantVars=relevantVars+1;
end
plotVars=index(1:relevantVars);


if printFlag==1
    figure
    for i=plotVars'
        scatter(A(:,i,N_runs),yA(:,N_runs))
        title(['x',num2str(i),'-y']);
        xlabel(['x',num2str(i)]);
        ylabel('y');
        figfile = fullfile(pathSAsub, ['X',num2str(i),'_Y_scatters.png']);
        saveas(gcf,figfile)
        
        
        
        
        
        figure
        histogram(A(:,i,N_runs))
        title(['x',num2str(i)])
        figfile = fullfile(pathSAsub, ['X',num2str(i),'_histogram.png']);
        saveas(gcf,figfile)
        
    end
    
    figure
    histogram(yA(:,N_runs))
    title('y')
    figfile = fullfile(pathSAsub, 'response_histogram.png');
    saveas(gcf,figfile)
    
    [f,x] = ecdf(yA(:,N_runs));
    figure
    plot(x,f)
    xlabel('v');
    ylabel('P(y<v)');
    figfile = fullfile(pathSAsub, 'response_ecdf.png');
    saveas(gcf,figfile)
    
    
    
if ~(size(samples,2)==1)
    figure
    for i=plotVars'
        hold on
        semilogx(samples,ST(i,1:N_runs))
    end
    xlabel('Sample Size');
    ylabel('ST');
    figfile = fullfile(pathSAsub, 'ST_convergence.png');
    saveas(gcf,figfile)
end


end


toc

% for i=1:3
% figure
% hold on
% plot(samples,SV(i,:))
% plot(samples,repmat(analyticalS(i), size(samples)))
% end
