function [ x, fval, exitflag, output ] = optimise_fmincon_wrapper(problem,par)

n_agents = 1;
if (isfield(par,'initial_population') && ~isempty(par.initial_population))
    initial_population = par.initial_population;
else
    if(isfield(par,'n_agents') && ~isempty(par.n_agents))
    % Initialise populations
        n_agents = par.n_agents;
    else
        n_agents = 1;
    end
    
    initial_population = lhsgen(n_agents,problem.dim).*repmat(problem.ub-problem.lb,n_agents,1)+repmat(problem.lb,n_agents,1);
    
end


% Run
f_best = nan;
x = nan(1,problem.dim);
fval = nan;
exitflag = nan;
output = struct;
i_best = 0;
options = optimset('Display','none',...%'MaxFunEvals',50*problem_fix_d.dim,'TolFun',1e-8,...%'LargeScale','off',...
                    'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
for i=1:size(initial_population,1)
    [x_i,f_i,exitflag_i,output_i] = fmincon(problem.objfun,initial_population(i,:),[],[],[],[],problem.lb,problem.ub,[],options,problem.par_objfun); %unconstrained
    if (i==1 || f_i < f_best)
        x = x_i;
        fval = f_i;
        exitflag = exitflag_i;
        output = output_i;
        i_best = i;

        fbest = f_i;
    end
end
output.i_best = i_best;

return