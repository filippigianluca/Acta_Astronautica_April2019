% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,par)
%
% input parameters:
%                  -) problem.dim
%                  -) problem.lb
%                  -) problem.ub
%                  -) problem.fitnessfcn
%
%                  -) par.n_agents
%                  -) par.n_populations





if (isfield(par,'initial_population') && ~isempty(par.initial_population))
    initial_population = par.initial_population;
else
    % Initialise populations
    initial_population = zeros(par.n_agents,problem.dim,par.n_populations);

    for s = 1 : par.n_populations
        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
        pop = lhsgen(par.n_agents,problem.dim).*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);

        
        initial_population(:,:,s) = pop;
    end
end


fitnessfcn.obj     = problem.objfun;
fitnessfcn.constr  = problem.constraint;
%--------------------------------------------------------------------------
% default parameters for MPAIDEA
fitnessfcn.obj_constr = par.fitnessfcn.obj_constr; %1;     % Flag to 1 if objective and constraints are in the same function
fitnessfcn.weighted   = par.fitnessfcn.weighted;   %0;     % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
fitnessfcn.ceq_eps    = par.fitnessfcn.ceq_eps;    %1e-6;  % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
fitnessfcn.w_ceq      = par.fitnessfcn.w_ceq;      %1000;  % Weights for penalty
fitnessfcn.w_c        = par.fitnessfcn.w_c;        %100;
%--------------------------------------------------------------------------




% Run MP-AIDEA
[memories_record, memories, archivebest, archiveALL, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, problem.lb, problem.ub, initial_population, par, problem.par_objfun);



% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), problem.dim);
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
%     x(i,:)  = memories_record(end, 1:end-1, i);
%     fval(i) = memories_record(end, end, i);
    x(i,:)  = memories_record(end, 1:end-2, i);
    fval(i) = memories_record(end, end-1, i);

end



output.memories_record      = memories_record;
output.memories             = memories;
output.archivebest          = archivebest;
output.archiveALL           = archiveALL;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;
output.nfeval               = par.nFeValMax;
return