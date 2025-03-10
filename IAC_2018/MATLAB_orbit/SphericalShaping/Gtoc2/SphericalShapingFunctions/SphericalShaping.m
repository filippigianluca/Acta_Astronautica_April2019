function [thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime1,Ucart,Vcart] = SphericalShaping(DepBody,ArrBody,t_dep,TOF,nr,toll,thetaStep,maxIter,a2,h,plotIter,Vinf,alpha,beta)
%25/09/15
%Function that implements the Spherical Shaping method (Novak, Vasile-2011)
% MODIFIED VERSION FOR GTOC (Adds an excess velocity ad departure)
%------------------------------------------------------------------------%
%INPUT
%DepBody, ArrBody: structures with data on departure and arrival planets
%t_dep: initial time (MJD)
%TOF: time of flight (days)
%nr: number of revolutions
%toll: relative tolerance on the cycle for the TOF
%thetaStep: step of the angle theta for the computation of the integral
%maxIter: maximum number of iterations of the Newton loop
%a2: coefficient a2 of the spherical shaping
%h: step for the finite derivative approximation
%plotIter: flag to indicate if to plot or not the convergence history
%Vinf, alpha, beta = initial excess velocity (modulus, azimuth, elevation
%respect to t-n-h frame; velocity in AU/day)

%OUTPUT
%thetaVec: vector containing angle theta from departure to arrival
%R: radius vector of trajectory
%Phi: vector with angle Phi of trajectory
%u_mat: matrix containing the control at each moment of time
%DV: DeltaV of the trajectory [km/s]
%isLoopConverged: indicator of convergence of the loop
%cont: contator of the number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu_S
%% Obtain Departure and Arrival position and velocities of spacecraft, that correspond to Departure Body and Arrival Body
% Transforn the initial and arrival Kepler elements in spherical
% coordinates

% Find position and velocity of Departure Body
M_dep = DepBody.M0 + DepBody.n*(t_dep-DepBody.t0);
Theta_dep = M2theta(M_dep,DepBody.e);

% Transform the orbital elements into cartesian coordinates
[dep_rC,dep_vC] = KeplElem2rv(DepBody.a,DepBody.e,DepBody.i,DepBody.w,DepBody.Omega,Theta_dep,mu_S);

% Calculate the excess velocity in cartesian reference frame
if Vinf == 0
    VinfCart = 0;
else
    [VinfCart] = addExcessVelocity(Vinf,alpha,beta,DepBody,Theta_dep);
end

% Transform the cartesian coordinates into spherical coordinates (from
% Roegiers Thesis)
[dep_rS,dep_vS] = Cartesian2Spherical(dep_rC, dep_vC + VinfCart);

%------------------------------------------------------------------------%
% Find position and velocity of Arrival Body
M_arr = ArrBody.M0 + ArrBody.n*(t_dep+TOF-ArrBody.t0);
Theta_arr = M2theta(M_arr,ArrBody.e);

% Transform the orbital elements into cartesian coordinates
[arr_rC,arr_vC] = KeplElem2rv(ArrBody.a,ArrBody.e,ArrBody.i,ArrBody.w,ArrBody.Omega,Theta_arr,mu_S);

% Transform the cartesian coordinates into spherical coordinates (from
% Roegiers Thesis)
[arr_rS,arr_vS] = Cartesian2Spherical(arr_rC,arr_vC);

%------------------------------------------------------------------------%
% Check the conditions on the angles
if (dep_rS(3) < -pi/2) || (dep_rS(3) > pi/2)
    disp('ERROR - Phi_dep is out of range');
end

if (arr_rS(3) < -pi/2) || (arr_rS(3) > pi/2)
    disp('ERROR - Phi_arr is out of range');
end

% Check that theta_arrival is not smaller than theta_dep, to avoid the
% design of a retrograde orbit
if arr_rS(2) < dep_rS(2)
    arr_rS(2) = arr_rS(2) + 2*pi;
end

%------------------------------------------------------------------------%
% Write the boundary conditions as in (Novak, Vasile - 2011)

%Initial
R_i = dep_rS(1);
Theta_i = dep_rS(2);
Phi_i = dep_rS(3);
TPrime_i = (R_i*cos(Phi_i))/dep_vS(2);   % (1/ThetaPrime_i)
RPrime_i = dep_vS(1)*TPrime_i;
PhiPrime_i = dep_vS(3)*TPrime_i/R_i;

%Final
R_f = arr_rS(1);
Theta_f = arr_rS(2) + 2*pi*nr;
Phi_f = arr_rS(3);
TPrime_f = (R_f*cos(Phi_f))/arr_vS(2);   % (1/ThetaPrime_f)
RPrime_f = arr_vS(1)*TPrime_f;
PhiPrime_f = arr_vS(3)*TPrime_f/R_f;

%Calculate alpha_i/f
alpha_i = - (RPrime_i*PhiPrime_i)/(PhiPrime_i^2 + cos(Phi_i)^2);
alpha_f = - (RPrime_f*PhiPrime_f)/(PhiPrime_f^2 + cos(Phi_f)^2);

%Calculate C_i/f
C_i = -(mu_S*TPrime_i^2)/R_i^2 + 2*RPrime_i^2/R_i + R_i*(PhiPrime_i^2 + cos(Phi_i)^2) + ...
      - RPrime_i*PhiPrime_i*(sin(Phi_i)*cos(Phi_i))/(PhiPrime_i^2 + cos(Phi_i)^2);
C_f = -(mu_S*TPrime_f^2)/R_f^2 + 2*RPrime_f^2/R_f + R_f*(PhiPrime_f^2 + cos(Phi_f)^2) + ...
      - RPrime_f*PhiPrime_f*(sin(Phi_f)*cos(Phi_f))/(PhiPrime_f^2 + cos(Phi_f)^2);
%% Start the cycle to satisfy the boundary conditions for the transfer

% Define the vector thetaVec that defines the nodes in which the functions
% values are calculated
% The function linspace has to be used, otherwise the final element of the
% vector will be different from theta_final
% The step of the vector is not thetaStep, but thetaStep/2; this is needed
% for the numerical quadrature to compute DV and u
% Moreover better to have an odd number of elements in the vector
nElem = ceil( (Theta_f - Theta_i)/ (thetaStep/2) ) + 1;  %should go a -1
if (mod(nElem,2) == 0)
    nElem = nElem + 1;
end
thetaVec = linspace(Theta_i,Theta_f,nElem);
trueThetaStep = thetaVec(3)-thetaVec(1);        % With the notations of quadrature formulas, this is h, not h/2!!!

% Set the vectors of coefficients x = [a0,a1,a3,..,a6,b0,..,b3]
x = zeros(10,1);

% Create the vector B (which depends only upon the initial conditions
B = [1/R_i; 1/R_f; Phi_i; Phi_f; -RPrime_i/R_i^2; -RPrime_f/R_f^2; PhiPrime_i; PhiPrime_f; C_i - 2*RPrime_i^2/R_i; C_f - 2*RPrime_f^2/R_f];

% Construct the matrix A with the coefficients of the parameters (given
% from the evaluation of the base functions in the initial and final
% points)
    A = [1 Theta_i cos(Theta_i) Theta_i*cos(Theta_i) sin(Theta_i) Theta_i*sin(Theta_i) 0 0 0 0;
        1 Theta_f cos(Theta_f) Theta_f*cos(Theta_f) sin(Theta_f) Theta_f*sin(Theta_f) 0 0 0 0;
        0  0  0  0  0  0  cos(Theta_i) Theta_i*cos(Theta_i) sin(Theta_i) Theta_i*sin(Theta_i);
        0  0  0  0  0  0  cos(Theta_f) Theta_f*cos(Theta_f) sin(Theta_f) Theta_f*sin(Theta_f);
        0  1  -sin(Theta_i) (cos(Theta_i)-Theta_i*sin(Theta_i)) cos(Theta_i) (Theta_i*cos(Theta_i)+sin(Theta_i)) 0 0 0 0;
        0  1  -sin(Theta_f) (cos(Theta_f)-Theta_f*sin(Theta_f)) cos(Theta_f) (Theta_f*cos(Theta_f)+sin(Theta_f)) 0 0 0 0;
        0  0  0  0  0  0  -sin(Theta_i) (cos(Theta_i)-Theta_i*sin(Theta_i)) cos(Theta_i) (Theta_i*cos(Theta_i)+sin(Theta_i));
        0  0  0  0  0  0  -sin(Theta_f) (cos(Theta_f)-Theta_f*sin(Theta_f)) cos(Theta_f) (Theta_f*cos(Theta_f)+sin(Theta_f));
        0  0  R_i^2*cos(Theta_i) R_i^2*(2*sin(Theta_i)+Theta_i*cos(Theta_i)) R_i^2*sin(Theta_i) R_i^2*(Theta_i*sin(Theta_i)-2*cos(Theta_i))...
        -alpha_i*cos(Theta_i) -alpha_i*(2*sin(Theta_i)+Theta_i*cos(Theta_i)) -alpha_i*sin(Theta_i) -alpha_i*(Theta_i*sin(Theta_i)-2*cos(Theta_i));
        0  0  R_f^2*cos(Theta_f) R_f^2*(2*sin(Theta_f)+Theta_f*cos(Theta_f)) R_f^2*sin(Theta_f) R_f^2*(Theta_f*sin(Theta_f)-2*cos(Theta_f))...
        -alpha_f*cos(Theta_f) -alpha_f*(2*sin(Theta_f)+Theta_f*cos(Theta_f)) -alpha_f*sin(Theta_f) -alpha_f*(Theta_f*sin(Theta_f)-2*cos(Theta_f))];

%------------------------------------------------------------------------%
% Set other parameters and flags
err = toll+1;          %error on the time of flight difference
isDpos = 1;            %indicator for D>0
isLoopConverged = 0;   %indicator of loop convergence
cont = 0;              %contator for the loop

% Set the vectors for the plot of (a2,TOFcomp)
a20 = a2;
if plotIter==1
 A2v = [];
 Tcv = [];
end
TOFcomp = [];

% Set the parameters for the line search algorithm (to modify the step of
% the Newton method)
step = 1;
fixVal = 0.05;
ind = 0;

%------------------------------------------------------------------------%
% Start the cycle to obtain the coefficients
while (err > toll) && (cont < maxIter)
    
    % Construct the vector Aa2
    Aa2 = [a2*Theta_i^2; a2*Theta_f^2; 0; 0; a2*2*Theta_i; a2*2*Theta_f; 0; 0; -a2*2*R_i^2; -a2*2*R_f^2];
    
    % Solve for the parameters x
    x = A\(B-Aa2);
    
    % Calculate the functions and their derivatives
    [R,Phi,RPrime1,PhiPrime1,RPrime2,PhiPrime2,RPrime3,PhiPrime3,D,DPrime,isDpos] = baseFuncDeriv(thetaVec,[x;a2]);
    
    if isDpos==1
        % Calculate T' and T''
        [TPrime1,TPrime2] = timeFunc(R,D,RPrime1,DPrime);
        
        % Calculate the time of flight with a Simpson quadrature
        TOFcomp(cont+1) = (trueThetaStep/6)*(TPrime1(1)+2*sum(TPrime1(3:2:end-1))+4*sum(TPrime1(2:2:end-1))+TPrime1(end));
        err = abs(TOFcomp(cont+1) - TOF)/TOF;
    end
    
    % Calculate new parameter a2 (only if errr > toll) with Newton-Raphson
    % iterations (***CONSIDER TO MAKE THIS A LOOP TO REDUCE THE ERROR
    % BETWEEN TOFcomp AND TOF TO ZERO***)
    if (err > toll)  && (isDpos==1)
        
        DeltaTOF = TOFcomp(cont+1) - TOF;
        
        % Calculate derivative of DeltaTOF with a central difference scheme
        [Rp,~,RPrime1p,~,~,~,~,~,Dp,DPrimep,isDposp] = baseFuncDeriv(thetaVec,[x;a2+h]);
        [Rm,~,RPrime1m,~,~,~,~,~,Dm,DPrimem,isDposm] = baseFuncDeriv(thetaVec,[x;a2-h]);
        
        % Check on ifDposp and ifDposm to avoid to obtain an imaginary time
        % of flight (and consequently an imaginary a2)
        if isDposp==1
            [TPrime1p,~] = timeFunc(Rp,Dp,RPrime1p,DPrimep);
            Ta2p = (trueThetaStep/6)*(TPrime1p(1)+2*sum(TPrime1p(3:2:end-1))+4*sum(TPrime1p(2:2:end-1))+TPrime1p(end));
            
            if isDposm==1
                [TPrime1m,~] = timeFunc(Rm,Dm,RPrime1m,DPrimem);
                Ta2m = (trueThetaStep/6)*(TPrime1m(1)+2*sum(TPrime1m(3:2:end-1))+4*sum(TPrime1m(2:2:end-1))+TPrime1m(end));
                
                % Finite derivative approximation
                dDeltaTOF = (Ta2p-Ta2m)/(2*h);
                
                % Step size modification based on a fixed (guessed) value
                if cont > 20
                    if (TOFcomp(cont) ~= 0) && (TOFcomp(cont+1) ~= 0)
                        if TOFcomp(cont) > TOF
                            if  ( TOFcomp(cont+1) > TOFcomp(cont) )                                              
                                if dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        else
                            if  ( TOFcomp(cont+1) < TOFcomp(cont) )
                                if dDeltaTOF > 0                                           % If dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        end
                    end
                end
                
                a2 = a2 + step*DeltaTOF/dDeltaTOF;
                
            else
                dDeltaTOF = (Ta2p-TOFcomp(cont+1))/(h);
                
                % Step size modification based on a fixed (guessed) value
                % In case try to cancel this part
                if cont > 20
                    if (TOFcomp(cont) ~= 0) && (TOFcomp(cont+1) ~= 0)
                        if TOFcomp(cont) > TOF
                            if  ( TOFcomp(cont+1) > TOFcomp(cont) )                                              
                                if dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        else
                            if  ( TOFcomp(cont+1) < TOFcomp(cont) )
                                if dDeltaTOF > 0                                           % If dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        end
                    end
                end
                
                a2 = a2 + step*DeltaTOF/dDeltaTOF;
            end
        else
            if isDposm==1
                [TPrime1m,~] = timeFunc(Rm,Dm,RPrime1m,DPrimem);
                Ta2m = (trueThetaStep/6)*(TPrime1m(1)+2*sum(TPrime1m(3:2:end-1))+4*sum(TPrime1m(2:2:end-1))+TPrime1m(end));
                
                dDeltaTOF = (TOFcomp(cont+1) - Ta2m)/(h);
                
                % Step size modification based on a fixed (guessed) value
                % In case try to cancel this part
                if cont > 20
                    if (TOFcomp(cont) ~= 0) && (TOFcomp(cont+1) ~= 0)
                        if TOFcomp(cont) > TOF
                            if  ( TOFcomp(cont+1) > TOFcomp(cont) )                                              
                                if dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        else
                            if  ( TOFcomp(cont+1) < TOFcomp(cont) )
                                if dDeltaTOF > 0                                           % If dDeltaTOF > 0
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                    ind = 1;
                                else
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            elseif ( abs(TOFcomp(cont+1)-TOFcomp(cont))/TOFcomp(cont) < 0.01 )
                                if ind==0
                                    step = abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                else
                                    step = -abs(fixVal/(DeltaTOF/dDeltaTOF))*rand;
                                end
                            else
                                if ind==0
                                    step=1;
                                else
                                    step=-1;
                                end
                            end
                        end
                    end
                end
                
                a2 = a2 + step*DeltaTOF/dDeltaTOF;
            else
                % Not possible to approximate the derivative in this case
                isDpos = 0;
                a2 = (a20-h) + step*2*h*rand;
            end
        end
        
    elseif (err > toll) && (isDpos==0) % In this case can't use Newton loop to find next value of a2
        % because TOFcomp would be an
        % imaginary number (if D<0)
        % Obtain the next a2 with a random variation
        a2 = (a20-h) + step*2*h*rand;
        
        % Need to set also a TOFcomp, otherwise risk of crash
        TOFcomp(cont+1) = 0;
        
    elseif  (err <= toll) && (isDpos==1)
        
        % In this case the error is smaller than the tolerance and isDpos = 1
        % so the loop has converged
        isLoopConverged = 1;
        
    end
    
    % Update cont (will have effect only if new iterations are required)
    cont = cont + 1;
    
    % Save the coeff a2 and the TOF comp for plot
    if (isDpos==1) && (plotIter==1)
        A2v(cont) = a2;
        T2v(cont) = TOFcomp(cont);
    elseif (isDpos==0) && (plotIter==1)
        A2v(cont) = a20;
        T2v(cont) = TOF;
    end
end

% Out from the loop, now check results
%------------------------------------------------------------------------%
if isLoopConverged==1
    % In this case the loop has converged, is so possible to calculate the
    % parameters of interest
    
    % Calculate thetaDot1 and thetaDot2 (first and second derivative)
    thetaDot1 = 1./TPrime1;
    thetaDot2 = -TPrime2./(TPrime1.^3);
    
    % Compute the r,v,a vectors in the reference
    % frame radial,orthoradial,out of plane
    U = PhiPrime1.^2 + cos(Phi).^2;
    Vtilde_mat = [RPrime1;
        R.*sqrt(U);
        zeros(1,length(thetaVec))];    % eq.16 Novak_Vasile
    Atilde_mat = [RPrime2-R.*U;
        2*RPrime1.*sqrt(U)+R.*PhiPrime1.*(PhiPrime2-sin(Phi).*cos(Phi))./(sqrt(U));
        (cos(Phi).*(PhiPrime2-sin(Phi).*cos(Phi)) + 2*sin(Phi).*U).*R./(sqrt(U))];   % eq.17 Novak_Vasile
    
    % Calculate the matrices with the versors at each instant of time
    % (angle theta), for the reference frames r-o-h and t-n-h
    % For the unit vector, just define them, not calculate them as
    % eh = cross(er,v)/norm(cross(er,v)) and e0 = cross(eh,er), since V
    % and a are already expressed in the r-o-h frame
    
    Er = [1; 0; 0];
    Eo = [0; 1; 0];
    Eh = [0; 0; 1];
    
    Et_mat = [];
    En_mat = [];
    
    for i=1:length(thetaVec)
        Et_mat(:,i) = Vtilde_mat(:,i)/norm(Vtilde_mat(:,i));
        En_mat(:,i) = cross(Eh,Et_mat(:,i));
    end
    
    % Compute the control acceleration (in the t-n-h frame) (in AU/day^2)
    u_mat = [];
    uMod_vec = [];
    
    for i=1:length(thetaVec)
        u_mat(1,i) = (mu_S/R(i)^2)*dot(Er,Et_mat(:,i)) + thetaDot2(i)*dot(Vtilde_mat(:,i),Et_mat(:,i)) +...
            + thetaDot1(i)^2*dot(Atilde_mat(:,i),Et_mat(:,i));
        u_mat(2,i) = (mu_S/R(i)^2)*dot(Er,En_mat(:,i)) + thetaDot1(i)^2*dot(Atilde_mat(:,i),En_mat(:,i));
        u_mat(3,i) = thetaDot1(i)^2*dot(Atilde_mat(:,i),Eh);
        
        uMod_vec(i) = sqrt(u_mat(1,i)^2 + u_mat(2,i)^2 + u_mat(3,i)^2);
    end
    
    % Compute the DeltaV  (in AU/day)
    DV = (trueThetaStep/6)*(TPrime1(1)*uMod_vec(1) + 2*sum( TPrime1(3:2:end-1).*uMod_vec(3:2:end-1) ) + ...
        + 4*sum(TPrime1(2:2:end-1).*uMod_vec(2:2:end-1)) + TPrime1(end)*uMod_vec(end));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Quiver3 plot and optimal control verification
    % Part to obtain the acceleration vector from the r-o-h frame to the
    % inertial frame (passing for the spherical frame) in order to plot it
    % with quiver3

    Ucart = [];
    Vtilde_cart = [];
    Vcart = [];
    % Define rotation matrices from R to S and from S to C and finally
    % calculate U vector directly in cartesian frame
    for i=1:length(thetaVec)
        Uroh(1,i) = (mu_S/R(i)^2)*dot(Er,Er) + thetaDot2(i)*dot(Vtilde_mat(:,i),Er) + thetaDot1(i)^2*dot(Atilde_mat(:,i),Er);
        Uroh(2,i) = (mu_S/R(i)^2)*dot(Er,Eo) + thetaDot2(i)*dot(Vtilde_mat(:,i),Eo) + thetaDot1(i)^2*dot(Atilde_mat(:,i),Eo);
        Uroh(3,i) = (mu_S/R(i)^2)*dot(Er,Eh) + thetaDot2(i)*dot(Vtilde_mat(:,i),Eh) + thetaDot1(i)^2*dot(Atilde_mat(:,i),Eh);
        
        P_RS = [ 1               0                         0;
                 0     cos(Phi(i))/sqrt(U(i))    -PhiPrime1(i)/sqrt(U(i));
                 0      PhiPrime1(i)/sqrt(U(i))    cos(Phi(i))/sqrt(U(i))];     
        
        P_SC = [ cos(thetaVec(i))*cos(Phi(i))    -sin(thetaVec(i))   -cos(thetaVec(i))*sin(Phi(i));
                 sin(thetaVec(i))*cos(Phi(i))     cos(thetaVec(i))   -sin(thetaVec(i))*sin(Phi(i));
                       sin(Phi(i))                      0                   cos(Phi(i))          ];
             
        Ucart(:,i) = P_SC*P_RS*Uroh(:,i);  
        Vtilde_cart(:,i) = P_SC*P_RS*Vtilde_mat(:,i);  
        Vcart(:,i) = Vtilde_cart(:,i)*thetaDot1(i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif isLoopConverged==0
    R=0;
    Phi=0;
    u_mat=0;
    DV=0;
    TPrime1=0;
    Ucart=0;
    Vcart=0;
end

%------------------------------------------------------------------------%
% Plot the graph a2,TOFcomp respect to iterations
if plotIter==1
    figure;
    subplot(2,1,1);
    plot([1:cont],A2v,'b-d',1,a20,'co');
    ylabel('a2');
    subplot(2,1,2);
    plot([1:cont],T2v,'r-d',[1:cont],TOF*ones(1,cont),'g');
    ylabel('TOF comp');
    xlabel('Iterations');
    hold on;
    suptitle('Coeff a2 and TOF comp with iterations');
end

end

