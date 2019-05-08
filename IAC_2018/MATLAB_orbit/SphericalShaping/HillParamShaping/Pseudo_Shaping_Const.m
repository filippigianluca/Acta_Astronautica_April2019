function [c,ceq] = Pseudo_Shaping_Const(y,LVec,LStep,A,b,TOF,mu,Lin_Trig)
%
% Pseudo_Shaping_Const: constraints function for the Pseudo Shaping.
% The constraint is the TOF
%
% INPUT
% y: parameters from the optimizer (free parameters)
% LVec: vector with angle u of the trajectory
% lStep: step of angle u of trajectory
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
l3 = y(3);

%--------------------------------------------------------------------------
% Extract L0,Lf
%--------------------------------------------------------------------------
L_0 = LVec(1);
L_f = LVec(end);
dL = L_f - L_0;

if Lin_Trig
    %--------------------------------------------------------------------------
    % Build the vector c
    %--------------------------------------------------------------------------
    c = [0;
         l1*sin(dL);
         0;
         l2*sin(dL);
         l2;
         l2*cos(dL);
         0
         l3*sin(dL);
         l3;
         l3*cos(dL)];
    
else
    %--------------------------------------------------------------------------
    % Build the matrix A
    %--------------------------------------------------------------------------
    A = [1 1 0 0 0 0 0 0 0 0;
         1 exp(l1*dL) 0 0 0 0 0 0 0 0;
         0 0 1 1 0 0 0 0 0 0;
         0 0 1 exp(l2*dL) 0 0 0 0 0 0;
         0 0 0 0 1 1 0 0 0 0;
         0 0 0 0 1 exp(l2*dL) 0 0 0 0;
         0 0 0 0 0 0 1 1 0 0;
         0 0 0 0 0 0 1 exp(l3*dL) 0 0;
         0 0 0 0 0 0 0 0 1 1;
         0 0 0 0 0 0 0 0 1 exp(l3*dL)];
    
    c = 0;
end
 
%--------------------------------------------------------------------------
% Solve for x
%--------------------------------------------------------------------------
x = A\(b-c);

%--------------------------------------------------------------------------
% Find the matrix with the Pseudo elements
%--------------------------------------------------------------------------
Pseudo_Matrix = Pseudo_Shaped(y,LVec,x,Lin_Trig);

% Extract the elements
p = Pseudo_Matrix(:,1);
f = Pseudo_Matrix(:,2);
g = Pseudo_Matrix(:,3);
%h = Pseudo_Matrix(:,4);
%k = Pseudo_Matrix(:,5);

%--------------------------------------------------------------------------
% Calculate time derivative of L
%--------------------------------------------------------------------------
q = (1 + f.*cos(LVec) + g.*sin(LVec));
r = p./q;
LDot = sqrt(mu*p).*(1./r).^2;
TPrime = 1./LDot;

%--------------------------------------------------------------------------
% Calculate the time of flight with a Simpson quadrature
%--------------------------------------------------------------------------
TOFcomp = (LStep/6)*(TPrime(1) + 2*sum(TPrime(3:2:end-1)) + 4*sum(TPrime(2:2:end-1)) + TPrime(end));
err = abs(TOFcomp - TOF)/TOF;
  
%--------------------------------------------------------------------------
% Set constraints
%--------------------------------------------------------------------------
c = [];
ceq = err;

end

