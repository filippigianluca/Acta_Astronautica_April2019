% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [ x, fval, exitflag, output ] = optimise_fmincon_wrapper(problem,par)
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
    x0 = par.initial_population;
else
    % Initialise populations

        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
        x0 = lhsgen(1,problem.dim).*repmat(problem.ub-problem.lb,1,1)+repmat(problem.lb,1,1);
    
end



% Run fmincon
[x,fval,eflag,outpt] = fmincon(problem.objfun, x0,[],[],[],[], problem.lb, problem.ub, problem.constraint, [], problem.par_objfun);



exitflag      = eflag;
output.nfeval = outpt.funcCount;

return