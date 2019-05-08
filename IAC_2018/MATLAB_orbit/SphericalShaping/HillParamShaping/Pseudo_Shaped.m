function Pseudo_Elements = Pseudo_Shaped(lambda,LVec,x,Lin_Trig)
%
% Pseudo_Shaped: function that calculates the Pseudo parameters correspondent
% to the boundary conditions,lambda parameters and independent variable l
%
% INPUT
% lambda: parameters from the optimizer (free parameters)
% LVec: vector with angle l of the trajectory
% x: vector with the coefficients of the boudnary conditions
% Lin_Trig: flag to select if to use Linear Trigonometric form or
% exponential form
%
% OUTPUT
% Pseudo_Elements: matrix with the Hill elements (n(uVec) x 6)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the matrix with the Pseudo elements

%--------------------------------------------------------------------------
% Extract from X the parameters _0, _1
%--------------------------------------------------------------------------
p0 = x(1);
p1 = x(2);
f0 = x(3);
f1 = x(4);
g0 = x(5);
g1 = x(6);
h0 = x(7);
h1 = x(8);
k0 = x(9);
k1 = x(10);

% Extract L_0
L_0 = LVec(1);

% Convert LVec in column vector
LVec = LVec(:);

%--------------------------------------------------------------------------
% Compute the Hill elements and store them in a matrix
%--------------------------------------------------------------------------

if Lin_Trig
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LINEAR TRIGONOMETRIC FORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %p
    p = p0 + p1*(LVec - L_0) + lambda(1)*sin(LVec - L_0);
    
    %f
    f = f0 + f1*(LVec - L_0) + lambda(2)*sin(LVec - L_0);
    
    %g
    g = g0 + g1*(LVec - L_0) + lambda(2)*cos(LVec - L_0);
    
    %h
    h = h0 + h1*(LVec - L_0) + lambda(3)*sin(LVec - L_0);
    
    %k
    k = k0 + k1*(LVec - L_0) + lambda(3)*cos(LVec - L_0);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXPONENTIAL FORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %p
    p = p0 + p1*exp(lambda(1)*(LVec - L_0));
    
    %f
    f = f0 + f1*exp(lambda(2)*(LVec - L_0));
   
    %g
    g = g0 + g1*exp(lambda(2)*(LVec - L_0));
    
    %h
    h = h0 + h1*exp(lambda(3)*(LVec - L_0));
    
    %k
    k = k0 + k1*exp(lambda(3)*(LVec - L_0));
end

%--------------------------------------------------------------------------
% Storage in Matrix
%--------------------------------------------------------------------------
Pseudo_Elements = [p f g h k];
end

