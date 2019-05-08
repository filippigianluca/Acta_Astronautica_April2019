function Hill_Elements = Hill_Shaped(lambda,uVec,x,Lin_Trig)
%
% Hill_Shaped: function that calculates the Hill parameters correspondent
% to the boundary conditions,lambda parameters and independent variable u
%
% INPUT
% lambda: parameters from the optimizer (free parameters)
% uVec: vector with angle u of the trajectory
% x: vector with the coefficients of the boudnary conditions
% Lin_Trig: flag to select if to use Linear Trigonometric form or
% exponential form
%
% OUTPUT
% Hill_Elements: matrix with the Hill elements (n(uVec) x 6)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the matrix with the Hill elements

%--------------------------------------------------------------------------
% Extract from X the parameters _0, _1
%--------------------------------------------------------------------------
a0 = x(1);
a1 = x(2);
a3 = x(3);
a4 = x(4);
G0 = x(5);
G1 = x(6);
H0 = x(7);
H1 = x(8);
h0 = x(9);
h1 = x(10);

% Extract u_0
u_0 = uVec(1);

% Convert uVec in column vector
uVec = uVec(:);

%--------------------------------------------------------------------------
% Compute the Hill elements and store them in a matrix
%--------------------------------------------------------------------------

if Lin_Trig
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LINEAR TRIGONOMETRIC FORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %r
    r = 1./(a0 + a1*uVec + lambda(3)*uVec.^2 + a3*cos(uVec) + a4*sin(uVec));
    
    %p
    p = -(r.^2).*(a1 + 2*lambda(3)*uVec - a3*sin(uVec) + a4*cos(uVec));
    
    %G
    G = G0 + G1*(uVec - u_0) + lambda(1)*sin(uVec - u_0);
    
    %H
    H = H0 + H1*(uVec - u_0) + lambda(2)*cos(uVec - u_0);
    
    %h
    h = h0 + h1*(uVec - u_0) + lambda(4)*sin(uVec - u_0);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXPONENTIAL FORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %r
    r = 1./(a0 + a1*uVec + lambda(3)*uVec.^2 + a3*cos(uVec) + a4*sin(uVec));
    
    %p
    p = -(r.^2).*(a1 + 2*lambda(3)*uVec - a3*sin(uVec) + a4*cos(uVec));
    
    %G
    G = G0 + G1*exp(lambda(1)*(uVec - u_0));
    
    %H
    H = H0 + H1*exp(lambda(1)*(uVec - u_0));
    
    %h
    h = h0 + h1*exp(lambda(2)*(uVec - u_0));
end

%--------------------------------------------------------------------------
% Storage in Matrix
%--------------------------------------------------------------------------
Hill_Elements = [r p G H h];
end

