% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [positions] = var2opt(i,problem_minmax)

dim_u_i = problem_minmax.dim_u_i;

positions_sorted = sum(dim_u_i(1:i-1))+1 : sum(dim_u_i(1:i)); 

positions = problem_minmax.order_dim_u(positions_sorted);
end