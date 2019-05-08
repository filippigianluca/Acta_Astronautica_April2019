function [uVec,timeVec,Hill_Matrix,R,V,ACont,DV,isLoopConverged] = Hill_Shaping_FMINCON(mu,Hill_Dep,Hill_Arr,t_dep,TOF,nr,toll,uStep,lambda0,Lin_Trig)
%
% Hill_Shaping_FMINCON: Function that realizes the shaping of a low thrust trajectory with the
% Hill parameters
%
% INPUT
% mu: gravitational parameter
% Hill_Dep: departure Hill elements
% Hill_Arr: arrival Hill elements
% t_dep: departure time
% TOF: time of flight
% nr: number of revolutions
% toll: tolerance
% uStep: step of angle u of trajectory
% lambda0: initial guess for the lambda (the optimizing variables)
% Lin_Trig: flag to select if to use Linear Trigonometric form or
% exponential form

% OUTPUT
% uVec: vector with angle u
% timeVec: vector with time
% Hill_Matrix: matrix with the Hill elements
% R,V: matrices with R,V
% ACont: matrix with control
% DV: deltaV value
% isLoopConverged: indicator on loop convergence
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Departure and Arrival Positions

%--------------------------------------------------------------------------
% Check on u
%--------------------------------------------------------------------------
u_0 = Hill_Dep(4);
u_f = Hill_Arr(4);

if (u_f < u_0)
    u_f = u_f + 2*pi;
end

% Add the number of revolutions
u_f = u_f + 2*pi*nr;

%% Define the Shaping parameters and variables

%--------------------------------------------------------------------------
% uVec: grid points for u (independet variable)
%--------------------------------------------------------------------------
% The step of the vector is not uStep, but uStep/2; this is needed
% for the numerical quadrature to compute DV and u
% Moreover better to have an odd number of elements in the vector
nElem = ceil( (u_f - u_0)/ (uStep/2) ) + 1;  %should go a -1
if (mod(nElem,2) == 0)
    nElem = nElem + 1;
end
uVec = linspace(u_0,u_f,nElem);
uVec = uVec(:);
trueUStep = uVec(3)-uVec(1);

%--------------------------------------------------------------------------
% Matrix A
%--------------------------------------------------------------------------
du = u_f - u_0;
A = [1 u_0 cos(u_0) sin(u_0) 0 0 0 0 0 0;
     1 u_f cos(u_f) sin(u_f) 0 0 0 0 0 0;
     0 1 -sin(u_0) cos(u_0) 0 0 0 0 0 0;
     0 1 -sin(u_f) cos(u_f) 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 1 du 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 1 du 0 0;
     0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 1 du];

%--------------------------------------------------------------------------
% Known vector b
%--------------------------------------------------------------------------
uDot_0 = Hill_Dep(3)/Hill_Dep(1)^2;
uDot_f = Hill_Arr(3)/Hill_Arr(1)^2; 

b = [1/Hill_Dep(1); 1/Hill_Arr(1);...
     -(Hill_Dep(2)/Hill_Dep(1)^2/uDot_0); -(Hill_Arr(2)/Hill_Arr(1)^2/uDot_f);...
     Hill_Dep(3); Hill_Arr(3);...
     Hill_Dep(5); Hill_Arr(5);...
     Hill_Dep(6); Hill_Arr(6)];

%--------------------------------------------------------------------------
% Solve with fmincon
%--------------------------------------------------------------------------
% Pass the initial guess to fmincon
optionsFM = optimoptions('fmincon','Algorithm','interior-point','Display','iter','TolX',1e-12,'TolCon',1e-8,'TolFun',1e-2,'MaxFunEvals',200000,'MaxIter',200);
[lambdaSol,fval,exitFlag] = fmincon(@(y) Hill_Shaping_Obj(y,uVec,trueUStep,A,b,mu,Lin_Trig,t_dep,TOF,toll),lambda0,[],[],[],[],[],[],@(y) Hill_Shaping_Const(y,uVec,trueUStep,A,b,TOF,Lin_Trig),optionsFM);

%--------------------------------------------------------------------------
% If solved, obtain the parameters
%--------------------------------------------------------------------------
if (exitFlag == 1) || (exitFlag ==2 )
    isLoopConverged = 1;
    
    % Compute solution
    l1 = lambdaSol(1);
    l2 = lambdaSol(2);
    a2 = lambdaSol(3);
    l4 = lambdaSol(4);
    
    if Lin_Trig
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
    
    x = A\(b-c);
    Hill_Matrix = Hill_Shaped(lambdaSol,uVec,x,Lin_Trig);
    
    r = Hill_Matrix(:,1);
    p = Hill_Matrix(:,2);
    G = Hill_Matrix(:,3);
    H = Hill_Matrix(:,4);
    h = Hill_Matrix(:,5);
    
    uDot = G./r.^2;
    TPrime = 1./uDot;

    TOFcomp = (uStep/6)*(TPrime(1) + 2*sum(TPrime(3:2:end-1)) + 4*sum(TPrime(2:2:end-1)) + TPrime(end));
    err = abs(TOFcomp - TOF)/TOF;
    fprintf('The TOFcomp is: %f and the real one is: %f with an error of %f\n',TOFcomp,TOF,err);
    
    DV = fval;
    
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
    % Ah = (G./r.^2 - uDot).*(G.*sin(i)./(r.*sin(uVec).*cos(i)));
    % %Ah = hDot.*(G.*sin(i)./(r.*sin(uVec)));
    % At = (HDot + G.*hDot.*(sin(i)).^2.*cot(uVec))./(r.*cos(i));
    %
    %
    % ACont(:,1) = Ar;
    % ACont(:,2) = At;
    % ACont(:,3) = Ah;
    %
    % % Modulus
    % Acc = sqrt( Ar.^2 + At.^2 + Ah.^2 );
      
else
    isLoopConverged = 0;
    Hill_Matrix = 0;
    ACont = 0;
    timeVec = 0;
    DV = 0;
    R = 0;
    V = 0;
    
end




end

