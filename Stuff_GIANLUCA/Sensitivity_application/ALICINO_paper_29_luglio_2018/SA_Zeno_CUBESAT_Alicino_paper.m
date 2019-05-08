clear
close all
clc


samples=[1e4];  


%--------------------------------------------------------------------------
objFun = @mask_SA_CUBESAT_ALICINO_paper_d_minmax;
[problem, algo_inner, algo_outer, algo_minmax, algo_decomposition] = def_CUBESAT_ALICINO();
for i=1:length(problem.lb_u{1})
    lb_u(i) = problem.lb_u{1}{i}(1);
    ub_u(i) = problem.ub_u{1}{i}(2);
end
lb_d = problem.lb_d;
ub_d = problem.ub_d;
%--------------------------------------------------------------------------


% N_x= length(lb_d) + length(lb_u);
% lb = [lb_d' lb_u];
% ub = [ub_d' ub_u];
par.dim_d = length(lb_d);


%----------------------
par.d = [10,90,7,0.390622153009282,0.0626265592897471,0.226873848743435,5,0.302065094536245,1];
N_x= length(lb_u);
lb = lb_u;
ub = ub_u;
%---------------------

%%%%%%%%%%%%%%%%%%%%%%
distributionType='Uniform';
order=1;
printFlag=1;
%%%%%%%%%%%%%%%%%%%%%%


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



figure
scatter(1:size(SV,1), SV(:,end),'filled')
title('SV')

figure
scatter(1:size(SV,1), ST(:,end),'filled')
title('ST')

[rank,index]=sort(ST(:,N_runs),'descend');










