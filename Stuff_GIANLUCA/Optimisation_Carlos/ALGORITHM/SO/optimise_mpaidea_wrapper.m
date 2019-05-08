function [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,par)

if (isfield(par,'initial_population') )
    if(size(par.initial_population,1) >= par.n_agents)
        initial_population = par.initial_population(1:par.n_agents,:,:);
    else
        initial_population = zeros(par.n_agents,problem.dim,par.n_populations);
        initial_population_0 = par.initial_population;
        remaining = par.n_agents - size(initial_population_0,1);
        for s = 1 : par.n_populations
        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
            pop = lhsgen(remaining,problem.dim).*repmat(problem.ub(1:problem.dim)-problem.lb(1:problem.dim),remaining,1)+repmat(problem.lb(1:problem.dim),remaining,1);
            initial_population(:,:,s) = [initial_population_0(:,:,s); pop];
        end
                
    end
else
    % Initialise populations
    initial_population = zeros(par.n_agents,problem.dim,par.n_populations);

    for s = 1 : par.n_populations
        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
        pop = lhsgen(par.n_agents,problem.dim).*repmat(problem.ub(1:problem.dim)-problem.lb(1:problem.dim),par.n_agents,1)+repmat(problem.lb(1:problem.dim),par.n_agents,1);

        
        initial_population(:,:,s) = pop;
    end
end

% Run MP-AIDEA
[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(problem.objfun, problem.lb(1:problem.dim), problem.ub(1:problem.dim), initial_population, par, problem.par_objfun);

% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), problem.dim);
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = memories_record(end, end, i);
end

output.memories             = memories;
output.archivebest          = archivebest;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;
output.nfeval               = options.nFeValMax;
output.memories_record      = memories_record;
return