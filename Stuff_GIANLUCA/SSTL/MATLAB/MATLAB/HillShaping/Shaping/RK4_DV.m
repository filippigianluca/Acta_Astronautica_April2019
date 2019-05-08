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

function [t,x] = RK4_DV (U,x0,F,r,G,i,F_n)

if size(x0,1)<size(x0,2) % Typical error, x0 given as row-vector
    x0=x0';
end

t = U ;
x = [x0] ;

lt=length(t);

for ii = 1:(lt-1)
    h = U(ii+1)-U(ii);
    k1 = DV(t(ii),x(ii),F,r,G,i,F_n);
    x1 = x(:,ii) + h/2*k1;
    k2 = DV(t(ii)+h/2,x1,F,r,G,i,F_n); 
    x2 = x(:,ii) + h/2*k2;
    k3 = DV(t(ii)+h/2,x2,F,r,G,i,F_n);
    x3 = x(:,ii) + h*k3;
    k4 = DV(t(ii+1),x3,F,r,G,i,F_n);
    x = [x , (x(:,ii) + h*(k1 + 2*k2 + 2*k3 + k4)/6) ];
    
end

end