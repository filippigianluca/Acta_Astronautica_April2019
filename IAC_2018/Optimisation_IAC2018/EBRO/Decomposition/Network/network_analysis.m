function problem = network_analysis(problem, minmax, minmin)
%
% input:  problem * CM  -> Coupled Matrix of uncertain vecor
%
% output: problem * CIM -> Coupled Indicator Matrix: matrix of the position
%                          with respect to the input u-vector of u-coupled
%                          and u-uncoupled

CM       = problem.CM;
[nr, nc] = size(CM);

n_obj    = problem.n_obj;
position_coupled_u = sum(problem.dim_u_i(1:problem.num_functions))+1:sum(problem.dim_u_i(:));


%% Coupled Indicator Matrix
CIM{nc,nc} = []; % matrix indicator of the coupled and uncoupled components
for i = 1 : nr
    nz = find(CM(i,:)==1);
    
    if length(nz)==1
        CIM{nz,nz} = [CIM{nz,nz}, i];
    elseif length(nz)==2
        CIM{nz(1),nz(2)} = [CIM{nz(1),nz(2)}, i];
    end
    
end



%% from CIM matrix to vector
dim_u_i       = [];
order_dim_u   = [];
matrix2linear = [];
U_Adj = zeros(nc, nc);

for j =1:nc
    order_dim_u = [order_dim_u CIM{j,j}];
    dim_u_i = [dim_u_i length(CIM{j,j})];
    matrix2linear = [matrix2linear {[j j]}];
end

for j =1:nc
    for k=1:nc
        if j<k
            %             if ~isempty(CIM{j,k})
            order_dim_u = [order_dim_u CIM{j,k}];
            dim_u_i = [dim_u_i length(CIM{j,k})];
            U_Adj(j,k) = 1;
            matrix2linear = [matrix2linear {[j k]}];
            %             end
        end
    end
end


%% Adjacency Matrix
Adj = U_Adj + U_Adj';



%% find nFE for u-vectors to sample
N_ind_components_u    = cell(1,1);
FE_u2sample           = cell(1,1);
FE_u2sampleVSuminmax  = [];
FE_u2sample_position  = [];
FE_minmax_uc          = [];

if ~isempty(problem.input.u_sample{n_obj})
    for i=1:size(problem.input.u_sample{1},1)
        FE_minmax_uc_k = [];
        FE_u2sample{i} = find_FE(problem.n_obj, problem.input.u_sample{1}(i,:), problem);
        FE_u2sample_position = [FE_u2sample_position; FE_u2sample{i}.position_FE];
        FE_u2sampleVSuminmax = [FE_u2sampleVSuminmax; minmax.FE.position_FE == FE_u2sample{i}.position_FE];
        
        for k = problem.num_functions +1 : length(problem.dim_u_i)
            FE_minmax_uc_k = [FE_minmax_uc_k min(FE_u2sampleVSuminmax(i,var2opt(k,problem)))];
        end
        
        FE_minmax_uc(i,:) = FE_minmax_uc_k;
    end
    
    for kkk = problem.num_functions+1 : length(problem.dim_u_i)
        [~, N_ind_components_u{kkk}] = unique(FE_u2sample_position(:,var2opt(kkk,problem)),'rows');
    end
end


%% matrix of the number of FEs for each coupling vector for all the sampling points
FE_u2Sample_M = cell(nc, nc);
AllSample = [minmax.u; problem.input.u_sample{1, 1}];
for j =1:nc
    for k=1:nc
        if k>j
            [~, b] = unique(AllSample(:,CIM{j,k}),'rows');
            FE_u2Sample_M(j,k) = {b};
        end
    end
end



%% plot graph proprierties
% dotMatrixPlot(Adj)
% drawCircGraph(Adj)


%% output
problem.Adj = Adj;
problem.CIM = CIM;

if isempty(problem.order_dim_u) && isempty(problem.dim_u_i)
    problem.order_dim_u = order_dim_u;
    problem.dim_u_i = dim_u_i;
end

for Nobj = 1:problem.n_obj
    for k=1:problem.dim_u
        problem.num_interval(k) = length(problem.bpa{Nobj}{k});
    end
end




problem.FE_sample_u            = FE_u2sampleVSuminmax;
problem.FE_minmax_uc           = FE_minmax_uc;
problem.N_ind_components_u     = N_ind_components_u;

problem.sample.FE_sample       = FE_u2sample;
problem.sample.FE_uSample_M    = FE_u2Sample_M;
problem.indexing.matrix2linear = matrix2linear;
return