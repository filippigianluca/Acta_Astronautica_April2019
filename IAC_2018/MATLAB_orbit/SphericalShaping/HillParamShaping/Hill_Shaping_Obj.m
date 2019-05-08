function Objective = Hill_Shaping_Obj(y,uVec,uStep,A,b,mu,Lin_Trig,t_dep,TOF,toll)
%
% Hill_Shaping_Obj: objective function for the Hill Shaping. Objective is
% DeltaV
%
%
% INPUT
% y: parameters from the optimizer (free parameters)
% uVec: vector with angle u of the trajectory
% uStep: step of angle u of trajectory
% A: matrix of the linear system to be solved
% b: known vector of boundary conditions of the linear system
% mu: gravitational parameter
% Lin_Trig: flag to select if to use Linear Trigonometric form or
% exponential form
% t_dep: departure time
% TOF: time of flight
% toll: tolerance

% OUTPUT
% objective: DeltaV
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the objective function (the DV)

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
p = Hill_Matrix(:,2);
G = Hill_Matrix(:,3);
H = Hill_Matrix(:,4);
h = Hill_Matrix(:,5);

%--------------------------------------------------------------------------
% Calculate time derivative of u
%--------------------------------------------------------------------------
uDot = G./r.^2;
TPrime = 1./uDot;

%--------------------------------------------------------------------------
% Compute position and velocity with transformation from Hill Elements to
% Kepler and then to Cartesian
%--------------------------------------------------------------------------
Kep_Matrix = Hill2Kep([Hill_Matrix(:,1:3) uVec Hill_Matrix(:,4:5)],mu);
Kep_Matrix = Kep_Matrix';

for i=1:length(uVec)
   [r,v] = KeplElem2rv(Kep_Matrix(i,1),Kep_Matrix(i,2),Kep_Matrix(i,3),...
       Kep_Matrix(i,5),Kep_Matrix(i,4),Kep_Matrix(i,6),mu);
   R(i,:) = r;
   %V(i,:) = v;
   RMod(i,1) = norm(r);
end

% Find the time evolution
[timeVec] = timeEvolution(uVec,TPrime,t_dep);
timeVec = timeVec';

% Velocity with finite approximation (time derivative)
for i=1:3
   V(:,i) = ( R(2:end,i) - R(1:end-1,i) )./( timeVec(2:end) - timeVec(1:end-1));
end

V(end+1,:) = V(end,:) + ( V(end,:) - V(end-1,:) );

    
% Acceleration with finite approximation (time derivative)
for i=1:3
   Acc(:,i) = ( V(2:end,i) - V(1:end-1,i) )./( timeVec(2:end) - timeVec(1:end-1));
end

Acc(end+1,:) = Acc(end,:) + ( Acc(end,:) - Acc(end-1,:) );

% Find the control acceleration
for i=1:3
   ACont(:,i) = Acc(:,i) + (mu./RMod.^3).* R(:,i);
end

% Compute modulus
AContMod = sqrt(ACont(:,1).^2 + ACont(:,2).^2 + ACont(:,3).^2);


% %--------------------------------------------------------------------------
% % Approximate the derivatives of the orbital elements with finite differences 
% %--------------------------------------------------------------------------
% pDot = ( p(2:end) - p(1:end-1) )./( uVec(2:end) - uVec(1:end-1) ).*uDot(1:end-1);
% pDot(end+1) = ( pDot(end) + pDot(end-1) )/2;
% 
% GDot = ( G(2:end) - G(1:end-1) )./( uVec(2:end) - uVec(1:end-1) ).*uDot(1:end-1);
% GDot(end+1) = ( GDot(end) + GDot(end-1) )/2;
% 
% HDot = ( H(2:end) - H(1:end-1) )./( uVec(2:end) - uVec(1:end-1) ).*uDot(1:end-1);
% HDot(end+1) = ( HDot(end) + HDot(end-1) )/2;
% 
% hDot = ( h(2:end) - h(1:end-1) )./( uVec(2:end) - uVec(1:end-1) ).*uDot(1:end-1);
% hDot(end+1) = ( hDot(end) + hDot(end-1) )/2;
% 
% %--------------------------------------------------------------------------
% % Calculate the control from the Gauss equations
% %--------------------------------------------------------------------------
% Ar = pDot - G.^2./r.^3 + mu./r.^2;
% An = (G./r.^2 - uDot).*(G.*sin(i)./(r.*sin(uVec).*cos(i)));
% %An = hDot.*(G.*sin(i)./(r.*sin(uVec)));
% At = (HDot + G.*hDot.*(sin(i)).^2.*cot(uVec))./(r.*cos(i));
% 
% % Modulus
% Acc = sqrt( Ar.^2 + An.^2 + At.^2 );


%--------------------------------------------------------------------------
% Calculate the DeltaV with a Simpson quadrature
%--------------------------------------------------------------------------
DV = (uStep/6)*(AContMod(1)*TPrime(1) + 2*sum(AContMod(3:2:end-1).*TPrime(3:2:end-1)) + ...
      + 4*sum(AContMod(2:2:end-1).*TPrime(2:2:end-1)) + AContMod(end)*TPrime(end));
  
%--------------------------------------------------------------------------
% Set Objective
%--------------------------------------------------------------------------
Objective = DV;

end

