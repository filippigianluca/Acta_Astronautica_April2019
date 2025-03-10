% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [minmax] = evaluate_minmax(problem_minmax, algo_minmax, algo_outer, algo_inner)






problem_minmax.sign_inner = 1;  %  for minmax



[ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax); % algo_minmax.optimise



minmax.d            = dmin;
minmax.u            = output.u;
minmax.f            = fminmax;
minmax.N_design_max = rank(dmin);
minmax.output       = output;

end