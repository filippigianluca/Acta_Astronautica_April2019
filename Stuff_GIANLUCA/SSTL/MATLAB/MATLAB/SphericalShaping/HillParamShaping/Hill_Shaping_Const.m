function [c,ceq] = Hill_Shaping_Const(y,uVec,uStep,A,b,TOF,Lin_Trig)
%
% Hill_Shaping_Const: constraints function for the Hill Shaping.
% The constraint is the TOF
%
% INPUT
% y: parameters from the optimizer (free parameters)
% uVec: vector with angle u of the trajectory
% uStep: step of angle u of trajectory
% A: matrix of the linear system to be solved
% b: known vector of boundary conditions of the linear system
% TOF: time of flight
% Lin_Trig: flag to select if to use Linear Trigonometric form or
% exponential form
%
% OUTPUT
% c: [] because no inequalities constraints
% ceq: equality constraint, on the TOF
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the constraint function (the TOF)

%--------------------------------------------------------------------------
% Extract lambda
%--------------------------------------------------------------------------
l1 = y(1);
l2 = y(2);
a2 = y(3);
l4 = y(4);

%--------------------------------------------------------------------------
% Extract u0,uf
%--------------------------------------------------------------------------
u_0 = uVec(1);
u_f = uVec(end);
du = u_f - u_0;

if Lin_Trig
    %--------------------------------------------------------------------------
    % Build the vector c
    %--------------------------------------------------------------------------
    c = [a2*u_0^2;
         a2*u_f^2;
         2*a2*u_0;
         2*a2*u_f;
         0;
         l1*sin(du);
         l2;
         l2*cos(du);
         0;
         l4*sin(du)];
    
else
    %--------------------------------------------------------------------------
    % Build the matrix A
    %--------------------------------------------------------------------------
    A = [1 u_0 cos(u_0) sin(u_0) 0 0 0 0 0 0;
        1 u_f cos(u_f) sin(u_f) 0 0 0 0 0 0;
        0 1 -sin(u_0) cos(u_0) 0 0 0 0 0 0;
        0 1 -sin(u_f) cos(u_f) 0 0 0 0 0 0;
        0 0 0 0 1 1 0 0 0 0;
        0 0 0 0 1 exp(l1*du) 0 0 0 0;
        0 0 0 0 0 0 1 1 0 0;
        0 0 0 0 0 0 1 exp(l1*du) 0 0;
        0 0 0 0 0 0 0 0 1 1;
        0 0 0 0 0 0 0 0 1 exp(l2*du)];
    
    c = [a2*u_0^2;
        a2*u_f^2;
        2*a2*u_0;
        2*a2*u_f;
        0;
        0;
        0;
        0;
        0;
        0];
   
end
 
%--------------------------------------------------------------------------
% Solve for x
%--------------------------------------------------------------------------
x = A\(b-c);

%--------------------------------------------------------------------------
% Find the matrix with the Hill elements
%--------------------------------------------------------------------------
Hill_Matrix = Hill_Shaped(y,uVec,x,Lin_Trig);

% Extract the elements
r = Hill_Matrix(:,1);
%p = Hill_Matrix(:,2);
G = Hill_Matrix(:,3);
%H = Hill_Matrix(:,4);
%h = Hill_Matrix(:,5);

%--------------------------------------------------------------------------
% Calculate time derivative of u
%--------------------------------------------------------------------------
uDot = G./r.^2;
TPrime = 1./uDot;

%--------------------------------------------------------------------------
% Calculate the time of flight with a Simpson quadrature
%--------------------------------------------------------------------------
TOFcomp = (uStep/6)*(TPrime(1) + 2*sum(TPrime(3:2:end-1)) + 4*sum(TPrime(2:2:end-1)) + TPrime(end));
err = abs(TOFcomp - TOF)/TOF;
  
%--------------------------------------------------------------------------
% Set constraints
%--------------------------------------------------------------------------
c = [];
ceq = err;

end

