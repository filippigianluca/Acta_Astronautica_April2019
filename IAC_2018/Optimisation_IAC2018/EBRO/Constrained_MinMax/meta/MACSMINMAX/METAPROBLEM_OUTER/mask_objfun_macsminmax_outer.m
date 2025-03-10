% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [varargout] = mask_objfun_macsminmax_outer(d,par_objfun)


[masked, ~, ~] = u_validation( par_objfun.problem_fix_d, d, par_objfun.u_record, par_objfun.local_search, par_objfun.objectives);


% if objfun and constraint are in the same function:
% varargout = [{f},{c},{ceq}]. Otherwise: varargout = f.
if isstruct(masked)  
    varargout{1} = masked.f;
    varargout{2} = masked.c;
    varargout{3} = masked.ceq;
else
    varargout{1} = masked;
end

return
