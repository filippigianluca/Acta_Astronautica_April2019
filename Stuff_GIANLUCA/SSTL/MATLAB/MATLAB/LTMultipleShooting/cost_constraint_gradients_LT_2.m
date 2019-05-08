function [myf,myc,myceq] =  cost_constraint_gradients_LT_2(x, parameters, options, constants)

% Function for the simultanous computation of objective function,
% constraints and their gradients
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

% Compute objective functions and constraints
[J,C,Ceq,GCeq] = cost_constraint_function_LT_2(x, parameters, options, constants);

myf.J = J;
myc.C = C;
myceq.Ceq = Ceq;

% If gradient of constraints or of objective functions is provided by the
% user, compute them with finite differences
if strcmp(options.options_fmincon.GradConstr, 'on') || ...
    strcmp(options.options_fmincon.GradObj, 'on')
    
    myceq.GCeq = GCeq';
  
else
    myf.GJ = [];
    myc.GC = [];
    myceq.GCeq = [];
end