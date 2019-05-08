% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [ x, fval, exitflag, output ] = optimise_constraint(problem,par)
% Maximisation of the constraint in the uncertain space.
% 
% min(-C)
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


% Run MP-AIDEA


% %% %%%%%%%%%%%%%%%%%%%
% problem.par_objfun.sign_constraint  = -1;
% %% %%%%%%%%%%%%%%%%%%%

% options = optimset('Display','none','MaxFunEvals',50*100,'TolFun',1e-8,...%'LargeScale','off',...
% 'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
% 
% [x,fval,~,output] = fmincon(@mask_constraints_macsminmax_max_constraint,initial_population,[],[],[],[],problem.lb, problem.ub, [],options, problem.par_objfun);




%--------------------------------------------------------------------------
% default parameters for MPAIDEA
fitnessfcn.obj_constr = 0;                         % Flag to 1 if objective and constraints are in the same function
fitnessfcn.weighted   = par.fitnessfcn.weighted;   %0;     % How to handle constraints: set to 1 for weighted constraints with fixed weights, or to 0 for penalty with no weights
fitnessfcn.ceq_eps    = par.fitnessfcn.ceq_eps;    %1e-6;  % If the constraints are handled without weights, then define a tolerance for the violation of the equality constraints
fitnessfcn.w_ceq      = par.fitnessfcn.w_ceq;      %1000;  % Weights for penalty
fitnessfcn.w_c        = par.fitnessfcn.w_c;        %100;
%--------------------------------------------------------------------------
fitnessfcn.obj       = @mask_constraints_macsminmax_max_constraint;  % Function to optimise
fitnessfcn.constr    = [];%problem.par_objfun.mask_constraints;    % Function of constraints
%--------------------------------------------------------------------------
problem.par_objfun.obj_constr = par.fitnessfcn.obj_constr;


[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, problem.lb, problem.ub, initial_population, par, problem.par_objfun);

%%%%%%%%%%%%%%%%%%%%%
% problem.par_objfun.sign_constraint  = 1;
%%%%%%%%%%%%%%%%%%%%%

% Output: minima and minima's objective value
x_all    = zeros(size(memories_record,3), problem.dim);
fval_all = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x_all(i,:)  = memories_record(end, 1:end-2, i);
    fval_all(i) = - memories_record(end, end -1, i);
end


[fval, n_constr] = max(fval_all);
x = x_all(n_constr, :);

output.memories             = memories;
output.archivebest          = archivebest;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;
output.nfeval               = par.nFeValMax;

return