function [x, f, eflag, outpt] = runobjconstr_LT(x0,parameters,options,constants,A,b,LB,UB)

% From Matlab help: Objective and Nonlinear Constraints in the Same Function
% 

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


fun = @objfun; % the objective function, nested below
cfun = @constr; % the constraint function, nested below

% Call fmincon
[x,f,eflag,outpt] = fmincon(fun,x0,A, b, [], [], LB, UB,cfun,options.options_fmincon);

    function [y,myGf] = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] =  cost_constraint_gradients_LT(x, parameters, options, constants);
        
            xLast = x;
        end
        
        y = myf.J;
        %         myGf = myf.GJ;
        myGf = zeros(length(x),1);
    end

    function [c,ceq,myGc,myGCeq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] =  cost_constraint_gradients_LT(x, parameters, options, constants);    
          
            xLast = x;
        end
        
        % Now compute constraint functions
        c = myc.C; % In this case, the computation is trivial
        ceq = myceq.Ceq;
        
        %         myGc = myc.GC;
        myGc = [];
        myGCeq = myceq.GCeq;
        
        
    end

end