function [t,x] = euler(ode_fun,t,x0,varargin)
%Euler integrator
%INPUT:  ode_fun = function of the ODE
%        t_limit = [t_in t_fin] vector with initial and final time
%        x0 = vector with initial integration point
%        h = step of the integration
%        varargin = in case more variables are needed to the ode_fun
%OUTPU:  t = vector with the evaluation times (all the points)
%        x = matrix [m*n], where m is the number of variables (the
%        dimension of x0, and n the number of time evaluation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation

%Calculate the length of the time vector (for the cycle)
N = length(t);

%Write the initial guess in the output matrix (as column vector)
x = x0(:);

%Write x_old for the cycle (and transpose)
x_old = x';

%Cycle for all the evaluation times
for i=2:N
    
    % Calculate the time step
    h = t(i) - t(i-1);
    
    %Calculate new vector x_new and add to the matrix x
    x_new = x_old + h*(ode_fun(t(i-1),x_old',varargin{:}))';
    x = [x x_new'];
    x_old = x_new;
    
end

end

