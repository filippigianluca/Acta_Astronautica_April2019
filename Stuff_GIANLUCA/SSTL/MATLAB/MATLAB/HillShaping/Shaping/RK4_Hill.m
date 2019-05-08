%% Cristian Greco 4520319 - RK4
% Input
% r_dot_control_function(t,x) = Function - Time derivative of the state x=r_dot_control_function(t,x);
% T = 2x1-vector - [T_initial T_final] Vector with initial and final 
% integration variable (usually time);
% h = Scalar - Step-size
% x0 = nx1-vector - Initial value of the state at T_initial

% Output - State integrated through RK-4
% x = matrix n x (1+(T_final-T_initial)/h) -  Column i = State at T_i
% t = (1+(T_final-T_initial)/h)-vector - Vector with epochs

function [t,x] = RK4_Hill (U,x0,F_r,F_u,F_n,mu)

if size(x0,1)<size(x0,2) % Typical error, x0 given as row-vector
    x0=x0';
end

t = U ;
x = [x0] ;

lt=length(t);

for ii = 1:(lt-1)
    h = U(ii+1)-U(ii);
    k1 = DiffHillEquationsWithControl(t(ii),x(:,ii),mu,F_r(t(ii)),F_u(t(ii)),F_n(t(ii)));
    x1 = x(:,ii) + h/2*k1;
    k2 = DiffHillEquationsWithControl(t(ii)+h/2,x1,mu,F_r(t(ii)+h/2),F_u(t(ii)+h/2),F_n(t(ii)+h/2)); 
    x2 = x(:,ii) + h/2*k2;
    k3 = DiffHillEquationsWithControl(t(ii)+h/2,x2,mu,F_r(t(ii)+h/2),F_u(t(ii)+h/2),F_n(t(ii)+h/2));
    x3 = x(:,ii) + h*k3;
    k4 = DiffHillEquationsWithControl(t(ii+1),x3,mu,F_r(t(ii+1)),F_u(t(ii+1)),F_n(t(ii+1)));
    x = [x , (x(:,ii) + h*(k1 + 2*k2 + 2*k3 + k4)/6) ];
    
end

end