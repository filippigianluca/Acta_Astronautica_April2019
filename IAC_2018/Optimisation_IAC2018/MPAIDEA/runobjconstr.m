function [x,f,eflag,outpt] = runobjconstr(x0, fc, flag_LG, opts, LB, UB, varargin)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% runobjconstr
% This function is used to:
% - compute objective function including weighted penalty for constrained
% problem to be handled by the DE in MPAIDEA
% - compute objective function for unconstrained problems for DE of MPAIDEA
% - compute solution of local search using fmincon, for both constrained
% and unconstrained problem
% Refer to: https://uk.mathworks.com/help/optim/ug/objective-and-nonlinear-constraints-in-the-same-function.html
%
%% Inputs:
%
% * x0: initial guess for fmincon, or individual of the population where to
%       evaluate the objective function
% * fc: structure containing information about the objective and
%       constraints, namely:
%       * fc.obj: objective function handle
%       * fc.constr: constraints function handle
%       * fc.w_ceq: weights for equality constraints
%       * fc.w_c: weights for inequality constraints
%       * fc.obj_constr: flag used to specify if the objectives and the
%                        constraints are defined in the same function
%                        (fc.obj_constr = 1) or in different functions
%                        (fc.obj_constr = 0)
%       * flag_LG: structure defining if this function should be used for
%                  the local search of MP-AIDEA (flag_LG.local = 1,
%                  flag_LG.global = 1) or for the DE evaluation in MP-AIDEA
%                  (flag_LG.local = 0, flag_LG.global = 1)
% * opts: options for optimisation with fmincon
% * LB   -> lower boundaries for fmincon
% * UB   -> upper boundaries for fmincon
% * varargin -> additional inputs

%% Outputs:
%
% * x: solution vector of fmincon if flag_LG.local=1 and flag_LG.global =
%      0; empty if flag_LG.local = 0 and flag_LG.global = 1
% * f: solution value of fmincon or weighted objective function for DE
% * eflag: exitflag for fmincon if flag_LG.local=1 and flag_LG.global = 0
% * outpt: output of fmincon if flag_LG.local=1 and flag_LG.global = 0

%% Author(s): Marilena Di Carlo (2015)
% email: marilena.di-carlo@strath.ac.uk


if nargin == 1 % No options supplied
    opts = [];
end

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


% -------------------------------------------------------------------------
% Local search: fmincon
% -------------------------------------------------------------------------
if flag_LG.local == 1
    
    fun = @objfun; % the objective function, nested below
    %     fun = @(xs)fc.obj(x,varargin);
    
    % Problem with upper and lower boundaries
    if ~isempty(LB) && ~isempty(UB)
        
        % If problem has constraints:
        if ~isempty(fc.constr) && ~fc.obj_constr
            cfun = @constr; % the constraint function, nested below
            % Call fmincon
            [x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],LB, UB,cfun,opts,varargin{:});
        elseif isempty(fc.constr) && fc.obj_constr
%             fun  = @obj_obj_constr;
%             cfun = @constr_obj_constr; % the constraint function, nested below
%             % Call fmincon
%             [x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],LB, UB,cfun,opts,varargin{:});           
            [x,f,eflag,outpt] = fmincon_runobjconstr(x0, opts, LB, UB, varargin{:});
        else
            % Call fmincon
            [x,f,eflag,outpt] = fmincon(fc.obj,x0,[],[],[],[],LB, UB,[],opts,varargin{:});
        end
        
        % Unbounded problem:
    else
        % Call fmincon
        [x,f,eflag,outpt] = fminunc(fun,x0,opts);
    end
    
    if ~isempty(fc.constr) || fc.obj_constr
        if eflag == -1 || eflag == -2 || outpt.constrviolation ~= 0
            
            if ~isempty(fc.constr) && ~fc.obj_constr
                constraint = feval(cfun,x0,varargin{:});
                func = feval(fun,x0,varargin{:});
            elseif isempty(fc.constr) && fc.obj_constr
                [func,constraint,~] = feval(fc.obj, x0,varargin{:});
            end
            
            if ( outpt.constrviolation ~= 0 && outpt.constrviolation > constraint )...
                    || ( outpt.constrviolation ~= 0 && outpt.constrviolation == constraint && func<f)  % feval(fun,x0,varargin{:})<f)
                x = x0;
                
                ff = feval(fun,x,varargin{:});% + outpt.constrviolation;
                
                if ~isempty(fc.constr) && ~fc.obj_constr
                    f = ff;
                elseif isempty(fc.constr) && fc.obj_constr
                    f = ff.f;
                end
                
                
                
            end
        end
    end
    % -------------------------------------------------------------------------
    % Global search DE: ojective function with penalty
    % -------------------------------------------------------------------------
elseif flag_LG.global == 1
    
    x = x0;
    
    % Evaluate objective function
    yyy = objfun(x, varargin{:});
    if ~isstruct(yyy)
        yy = yyy;
    end
    
    % If problem is constrained, evaluate function of constraints
    if ~isempty(fc.constr)
        [c,ceq] = constr(x, varargin{:});
    elseif isempty(fc.constr) && fc.obj_constr
        yy  = yyy.f;
        c   = yyy.c;
        ceq = yyy.ceq;
    else
        c = [];
        ceq = [];
    end
    
    if (fc.weighted && (~isempty(ceq) || ~isempty(c)) ) || ...
            ( isempty(ceq) && isempty(c) )
        
        % Defined weighted objective function:
        % Both equality and inequality constraints
        if ~isempty(ceq) && ~isempty(c)
            f = yy + fc.w_ceq * norm(ceq)  + fc.w_c * ( abs(c) .* (c>0 == 1));
            % Only equality constraints
        elseif ~isempty(ceq) && isempty(c)
            f = yy + fc.w_ceq * norm(ceq)  ;
            % Only inequality constraints
        elseif isempty(ceq) && ~isempty(c)
            f = yy + fc.w_c * ( abs(c) .* (c>0 == 1));
            % No constraints
        elseif isempty(ceq) && isempty(c)
            f = yy ;
        end
        
    elseif fc.weighted == 0 && (~isempty(ceq) || ~isempty(c))
        
        f.yy = yy;
        f.ceq = norm(ceq);
        f.c = c;
        
        if norm(f.ceq) > fc.ceq_eps  && any(f.c > 0)
            f.non_feas_ceq = 1;
            f.non_feas_c    = 1;
        elseif norm(f.ceq) > fc.ceq_eps  && all(f.c < 0)
            f.non_feas_ceq = 1;
            f.non_feas_c    = 0;
        elseif norm(f.ceq) < fc.ceq_eps  && any(f.c > 0)
            f.non_feas_ceq = 0;
            f.non_feas_c    = 1;
        else
            f.non_feas_ceq = 0;
            f.non_feas_c = 0;
        end
        
    end
    
    eflag = [];
    outpt = [];
    
end

    function y = objfun(x, varargin)
        
        % Check if computation is necessary (only if
        % objective function and constraint functions are the same)
        if ~isequal(x,xLast) && fc.obj_constr
            
            [myf,myc,myceq] = feval(fc.obj,x,varargin{:});
            xLast = x;
            y.f = myf;
            y.c = myc;
            y.ceq = myceq;
            
        elseif isequal(x,xLast) && fc.obj_constr
            y.f = myf;
            y.c = myc;
            y.ceq = myceq;
            
        elseif ~isequal(x,xLast) && ~fc.obj_constr
            % Evaluate objective function
            [myf] = feval(fc.obj,x,varargin{:});
            y = myf;
        end
        
        
    end

    function [c,ceq] = constr(x, varargin)
        % Check if computation is necessary
        if ~isequal(x,xLast) && fc.obj_constr
            
            [myf,myc,myceq] = feval(fc.constr, x,varargin{:});
            xLast = x;
            
        elseif isequal(x,xLast) && fc.obj_constr
            
        elseif ~isequal(x,xLast) && ~fc.obj_constr
            % Evaluate constraint function
            [myc,myceq] = feval(fc.constr, x,varargin{:});
            
        end
        
        % Now compute constraint functions
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end









%% new for constr and obj in the same function

% % only for fmincon when obj_constr == 1
%     function [c,ceq] = constr_obj_constr(x, varargin)
%         % Check if computation is necessary
%         [~,myc,myceq] = feval(fc.obj, x,varargin{:});
%         
%         
%         % Now compute constraint functions
%         c = myc;
%         ceq = myceq;
%     end
% 
% 
% % only for fmincon when obj_constr == 1
%     function y = obj_obj_constr(x, varargin)
%         % Check if computation is necessary
%         [myf] = feval(fc.obj, x,varargin{:});
%         
%         
%         % Now compute constraint functions
%         y = myf;
%     end



% only for fmincon when obj_constr == 1
function [x,f,eflag,outpt] = fmincon_runobjconstr(x0,opts,LB,UB,varargin)

if nargin == 1 % No options supplied
    opts = [];
end

fmincon_xLast = []; % Last place computeall was called
fmincon_myf   = []; % Use for objective at xLast
fmincon_myc   = []; % Use for nonlinear inequality constraint
fmincon_myceq = []; % Use for nonlinear equality constraint

fmincon_fun  = @fmincon_objfun;  % the objective function, nested below
fmincon_cfun = @fmincon_constr; % the constraint function, nested below


            
% Call fmincon
[x,f,eflag,outpt] = fmincon(fmincon_fun,x0,[],[],[],[],LB,UB,fmincon_cfun,opts,varargin{:});

    function y = fmincon_objfun(x,  varargin)
        if ~isequal(x,fmincon_xLast) % Check if computation is necessary
            [fmincon_myf,fmincon_myc,fmincon_myceq] = fc.obj(x, varargin{:});
            fmincon_xLast = x;
        end
%         % Now compute objective function
        y = fmincon_myf;
    end

    function [c,ceq] = fmincon_constr(x, varargin)
        if ~isequal(x,fmincon_xLast) % Check if computation is necessary
            [fmincon_myf,fmincon_myc,myceq] = fc.obj(x, varargin{:});
            fmincon_xLast = x;
        end
        % Now compute constraint functions
        c   = fmincon_myc; % In this case, the computation is trivial
        ceq = fmincon_myceq;
    end

end






end