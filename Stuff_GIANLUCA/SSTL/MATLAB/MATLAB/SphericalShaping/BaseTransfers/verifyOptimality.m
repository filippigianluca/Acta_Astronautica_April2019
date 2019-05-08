function verifyOptimality(DepBody,ArrBody,R,Vcart,Ucart,timeVec,mu,DV)
%
%verifyOptimality: function that verifies the optimality of a solution
%(X,Y,Z) with respect to the optimal control adjoint equations, for a
%problem of minimum DeltaV. Procedure: obtain the Lagrange multipliers in
%two way and compare them. Way 1: from the optimal control condition and
%with numerical derivative approximation. Way 2: from integration of
%adjoint differential equations, but starting with initial conditions given
%by the optimal control condition (because the true ones are not available).
%
% Very probably the result will not be accurate, because for calculating
% the adjoints with numerical integration, the time vector obtained from
% the spherical shaping is used; in that vector the time step is not
% regular and is quite big (on the order of days). The origin of this
% problem is the theta step given in input to the spherical shaping, which
% is used then to obtain all the subsequent vectors.
%
%In the last part the optimal trajectory is found using optimal control
%formulation, starting from initial conditions given by the spherical
%shaping method. Then a comparison between the two is done.
%
% The methods used are different: INDIRECT SINGLE SHOOTING, INDIRECT
% MULTIPLE SHOOTING, DIRECT TRANSCRIPTION, MATLAB bvp4c. There are also some variations,
% like scaled equations, linear control in the performance index, different
% discretization schemes and the
% use of the OPTI toolbox instead of Matlab functions. The selection of the
% method is done through a menu.
%
%
% INPUT
% DepBody, ArrBody: departure and arrival bodies
% R: radius vector in cartesian coordinates [x,y,z]. Is a matrix of size
% (3 x timeStep)
% Vcart: velocity in cartesian frame
% Ucart: control acceleration in cartesian frame
% timeVec: time vector of size (1 x timeStep)
% mu: gravitational parameter
%
% OUTPUT
% No output, just graph plotting
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Initialization

% Add path of subfolders containing the functions for optimality
%addpath(genpath('OptimalityFunctions'));
addpath('OptimalityFunctions');

%--------------------------------------------------------------------------
% CONSTANTS SETTING
%--------------------------------------------------------------------------
AU2km = 149597870.66;
km2AU = 1/AU2km;
T_sid = 86164.1004;

%--------------------------------------------------------------------------
% ALGORITHM SELECTION
%--------------------------------------------------------------------------
% Set flags to indicate which algorithm to perform
INDIRECT_SINGLE_SHOOTING = 0;
INDIRECT_MULTIPLE_SHOOTING = 0;
INDIRECT_MULTIPLE_SHOOTING_LIN_U = 0;
TPBVP = 0;
DIRECT_TRANSCRIPTION = 0;

% Select the Optimal Control Method
choice = menu('Select Optimal Control Method','INDIRECT_SINGLE_SHOOTING','INDIRECT_MULTIPLE_SHOOTING',...
    'INDIRECT_MULTIPLE_SHOOTING_LIN_U','TPBVP','DIRECT_TRANSCRIPTION','NO METHOD');
switch choice
    case 1
       INDIRECT_SINGLE_SHOOTING = 1;
    case 2
       INDIRECT_MULTIPLE_SHOOTING = 1;
    case 3
       INDIRECT_MULTIPLE_SHOOTING_LIN_U = 1; 
    case 4
       TPBVP = 1;
    case 5
       DIRECT_TRANSCRIPTION = 1; 
    case 6
end

%--------------------------------------------------------------------------
% FLAG SELECTION
%--------------------------------------------------------------------------
% Select (1) or not select (0) if to:
% 1) Integrate forward or backward the diff eq
% 2) Integrate only the adjoints diff eq or the complete diff eq
% 3) Plot the adjoints and Hamiltonian
% 4) Use the scaled version of the equations (where implemented)
% 5) Verify the control obtained by the selected direct or indirect method
% 6) Quadrature method for DV: 1 for CS composite, 0 for Gauss composite
% 7) Use the alternative OPTI Toolbox for solving (where implemented)
% 8) Save the trajectory obtained
% 9) Load an external trajectory
% 10) Compute a tangential thrust trajectory
FORWARD_INT = 0;
ONLY_ADJOINTS = 1;
PLOT_ADJOINTS_HAM = 1;
SCALED_EQ = 0;
VERIFY_CONTROL = 1;
CS_INT = 0;
OPTI = 0;
SAVE_TRAJ = 0;
LOAD_TRAJ = 0;
TANG_THRUST_TRAJ = 0;

%% 2) Calculate the history of the Langrange multipliers from the control
% obtained (multipliers correspondent to the shaped solution)

lambdavSS = -Ucart;
lambdavdotSS = [];

for i=1:3
    lambdavdotSS(i,:) = ( lambdavSS(i,2:end) - lambdavSS(i,1:end-1) )./( timeVec(2:end) - timeVec(1:end-1));
end

% Add last element (not available from finite differencing
% Use interpolation to extrapolate (doesnt work for now)
%lambdavdotSS(:,end+1) = interp1(timeVec(1:end-1),lambdavdotSS(:,:),timeVec(end),'nearest','extrap');
% Approximate as a mean of the previous two derivatives (very rough)
lambdavdotSS(:,end+1) = ( lambdavdotSS(:,end) + lambdavdotSS(:,end-1) )/2;

% Obtain lambda r
lambdarSS = -lambdavdotSS;

%% 3) Calculateh the history of the Lagrange multipliers from the adjoint differential equations

% Choose starting conditions depending upon integration direction
if (FORWARD_INT == 1)
    
    % Initial conditions for integration forward
    lambdar0 = lambdarSS(:,1);
    lambdav0 = lambdavSS(:,1);
    
    r0 = R(:,1);
    v0 = Vcart(:,1);
  
    % Create the Flip vector (even if in this case is not flipped)
    timeVecFlip = timeVec;
    
else
    
    % Initial conditions for integration backward
    lambdar0 = lambdarSS(:,end);
    lambdav0 = lambdavSS(:,end);
    
    r0 = R(:,end);
    v0 = Vcart(:,end);
    
    % Flip the time vector for backward integration
    timeVecFlip = fliplr(timeVec);
end

%--------------------------------------------------------------------------
% INTEGRATION OF ADJOINTS DIFFERENTIAL EQUATIONS
%--------------------------------------------------------------------------
% Integrate the differential equations (or only the adjoints or all the
% equations) forward or backward starting from the known boundary
% conditions (exact for r,v and estiamted for lambda). Recover the state at
% the boundary of integration and compare it to the previous state
% available (for r,v the exact state, for lambda the estimated)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use rk45 (Verified, also the backwards integration, with ode45)
%[t,lambdaAd] = rk4(@adjointEq, timeVec, [lambdar0;lambdav0], FORWARD_INT, R, timeVec, mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ONLY_ADJOINTS
    
    % Integrate only adjoints
    x0 = [lambdar0;lambdav0];
    [tOde,lambdaAdOde] = ode45(@(t,lambda) adjointEqOde(t,lambda,R,timeVec,mu), timeVecFlip, x0);
    
    % Reflip to obtain vectors ordered with increasing time
    if (FORWARD_INT == 1)
        t = tOde;
        lambdaAd = lambdaAdOde';
    else
        t = fliplr(tOde');
        lambdaAd = fliplr(lambdaAdOde');
    end
    
else
    
    % Integrate all the differential equations (only backward option) to see if
    % the initial state (r,v) is recovered correctly
    x0 = [r0;v0;lambdar0;lambdav0];
    [tOde,YOde] = ode45(@(t,x) diffEq(t,x,mu), timeVecFlip, x0);
    
    % Reflip to obtain vectors ordered with increasing time
    if (FORWARD_INT == 1)
        t = tOde;
        Y = YOde';
        lambdaAd = Y(7:12,:);
    else
        t = fliplr(tOde');
        Y = fliplr(YOde');
        lambdaAd = Y(7:12,:);
    end
    
    % Plot radius and velocity recovered from integration if the integration of
    % all the diff equations is performed
    
    if PLOT_ADJOINTS_HAM
        %Radius vector
        figure;
        subplot(3,1,1);
        plot(timeVec,R(1,:),'b','linewidth',2);
        hold on;
        plot(t,Y(1,:),'g','linewidth',2);
        ylabel('R x');
        legend('Inverse Method','Int From Adjoints','location','northeast');
        subplot(3,1,2);
        plot(timeVec,R(2,:),'b','linewidth',2);
        hold on;
        plot(t,Y(2,:),'g','linewidth',2);
        ylabel('R y');
        subplot(3,1,3);
        plot(timeVec,R(3,:),'b','linewidth',2);
        hold on;
        plot(t,Y(3,:),'g','linewidth',2);
        ylabel('R z');
        suptitle('Radius vector comparison Inverse Method-Integration From Adjoints');
        xlabel('Time');
        
        %Radius vector
        figure;
        subplot(3,1,1);
        plot(timeVec,Vcart(1,:),'b','linewidth',2);
        hold on;
        plot(t,Y(4,:),'g','linewidth',2);
        ylabel('V x');
        legend('Inverse Method','Int From Adjoints','location','northeast');
        subplot(3,1,2);
        plot(timeVec,Vcart(2,:),'b','linewidth',2);
        hold on;
        plot(t,Y(5,:),'g','linewidth',2);
        ylabel('V y');
        subplot(3,1,3);
        plot(timeVec,Vcart(3,:),'b','linewidth',2);
        hold on;
        plot(t,Y(6,:),'g','linewidth',2);
        ylabel('V z');
        suptitle('Velocity vector comparison Inverse Method-Integration From Adjoints');
        xlabel('Time');   
    end
    
end

%% 4) Compare the adjoints obtained from the optimal control condition and from the adjoint differential equations

%--------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------
if PLOT_ADJOINTS_HAM
    
    % R adjoints
    figure;
    subplot(3,1,1);
    plot(timeVec,lambdarSS(1,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(1,:),'g','linewidth',2);
    ylabel('\lambda rx');
    legend('Opt control condition','Adjoint diff equations','location','northeast');
    subplot(3,1,2);
    plot(timeVec,lambdarSS(2,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(2,:),'g','linewidth',2);
    ylabel('\lambda ry');
    subplot(3,1,3);
    plot(timeVec,lambdarSS(3,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(3,:),'g','linewidth',2);
    ylabel('\lambda rz');
    suptitle('r adjoints time evolution for Spherical Shaping');
    xlabel('Time');
    
    % V adjoints
    figure;
    subplot(3,1,1);
    plot(timeVec,lambdavSS(1,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(4,:),'g','linewidth',2);
    ylabel('\lambda vx');
    legend('Opt control condition','Adjoint diff equations','location','northeast');
    subplot(3,1,2);
    plot(timeVec,lambdavSS(2,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(5,:),'g','linewidth',2);
    ylabel('\lambda vy');
    subplot(3,1,3);
    plot(timeVec,lambdavSS(3,:),'b','linewidth',2);
    hold on;
    plot(t,lambdaAd(6,:),'g','linewidth',2);
    ylabel('\lambda vz');
    suptitle('v adjoints time evolution for Spherical Shaping');
    xlabel('Time');
    
    %--------------------------------------------------------------------------
    % HAMILTONIAN CALCULATION
    %--------------------------------------------------------------------------
    %Calculate the Hamiltonian at each moment of time
    % with the adjoints coming from the optimal control condition
    HSS = [];
    HAd = [];
    for i=1:length(timeVec)
        HSS(i) = 0.5*norm(Ucart(:,i))^2 + dot( lambdarSS(:,i),Vcart(:,i) ) + dot( lambdavSS(:,i), (  -mu*R(:,i)/norm(R(:,i))^3  + Ucart(:,i)) );
        HAd(i) = 0.5*norm(Ucart(:,i))^2 + dot( lambdaAd(1:3,i),Vcart(:,i) ) + dot( lambdaAd(4:6,i), (  -mu*R(:,i)/norm(R(:,i))^3  + Ucart(:,i)) );
    end
    
    % Hamiltonian plotting
    figure;
    plot(timeVec,HSS,'b','linewidth',2);
    hold on;
    plot(t,HAd,'g','linewidth',2);
    plot(timeVec,ones(1,length(timeVec))*HSS(1),'r','linewidth',2);
    xlabel('Time');
    ylabel('Hamiltonian function');
    suptitle('Hamiltionian with time');
    legend('H inv','H adj','Initial value');
    
end

%% Set BC and parameters

%--------------------------------------------------------------------------
% BOUNDARY CONDITIONS AND PARAMETERS
%--------------------------------------------------------------------------
% To be more correct, the initial and final positions and velocitites
% should be calculated from the Departure and Arrival bodies, but since 
% the inverse solution has been verified, it is known that such boundaries
% are exact (within tolerances)
t0 = timeVec(1);
tf = timeVec(end);
TOF = tf - t0;
r0 = R(:,1)';
v0 = Vcart(:,1)';
rf = R(:,end)';
vf = Vcart(:,end)';
uMax = DV/TOF;

%--------------------------------------------------------------------------
% LOADING OF REFERENCE TRAJECTORY (IF SELECTED)
%--------------------------------------------------------------------------
% Done after the trascription of boundary conditions because they
% must be exact (so coming from a validated solution, i.e. the inverse
% one), because there is no guarantee that a generic initial guess will
% respect the boundary conditions at the final time (tf)

% Selection
if LOAD_TRAJ
   dati = load('traj_rvut.mat'); 
   R = dati.r;
   Vcart = dati.v;
   Ucart = dati.Uopt;
   timeVec = dati.tVec;
end

%% 5) Calculate optimal trajectory with indirect single shooting

if INDIRECT_SINGLE_SHOOTING
    
    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  INDIRECT SINGLE SHOOTING %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % SELECTION OF LAMBDA
    %--------------------------------------------------------------------------
       lambda0(1:3) = lambdarSS(:,1)';
       lambda0(4:6) = lambdavSS(:,1)';
    %lambda0(1:3) = lambdaAd(1:3,1)';
    %lambda0(4:6) = lambdaAd(4:6,1)';
    %  lambda0(1:3) = lambdarSS(:,end)';
    %  lambda0(4:6) = lambdavSS(:,end)';
    %  lambda0(1:3) = lambdaAd(1:3,end)';
    %  lambda0(4:6) = lambdaAd(4:6,end)';
    
    %--------------------------------------------------------------------------
    % FUNCTION CALL
    %--------------------------------------------------------------------------
    % Call the function that computes the optimal trajectory
    % Results in AU,days
    [r,v,lambda,tVec] = OptimalTrajectoryIndirectFSOLVESingle(timeVec,t0,tf,r0,rf,v0,vf,lambda0,mu);
    
    %--------------------------------------------------------------------------
    % RESULTS PLOTTING
    %--------------------------------------------------------------------------
    
    % First obtain control from the optimal control condition
    Uopt = (-lambda(:,4:6))';
    
    % Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r,Uopt);
    
    % Plot comparison of r,v,u,lambda between inverse method and optimal
    % control
    comparisonLRVU(timeVec,tVec,R,r',Vcart,v',Ucart,Uopt,lambdaAd(1:3,:),lambda');

    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        % Transpose for compatibility
        r = r';
        v = v';
        save('traj_rvut.mat','r','v','Uopt','tVec');
        % Re-Transpose for compatibility
        r = r';
        v = v';
    end
    
    %--------------------------------------------------------------------------
    % DV CALCULATION
    %--------------------------------------------------------------------------
    if CS_INT
        % Use Cavalieri-Simpson composite quadrature rule (accurate only for
        % odd number of mesh points)
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        tStep = tVec(3) - tVec(1);
        DV = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
            + 4*sum(uMod(2:2:end-1)) + uMod(end));
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
        
    else  
        % Use Gauss 2 points composite quadrature rule (accurate for odd and
        % even number of mesh points).
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        alpha = (tVec(2)-tVec(1))/2;
        tVecN0 = ( tVec(1:end-1) + tVec(2:end) )/2 + alpha/(sqrt(3));
        tVecN1 = ( tVec(1:end-1) + tVec(2:end) )/2 - alpha/(sqrt(3));
        
        % Find uMod1 and uMod2 with an interpolation (is not possible, or too
        % complicated, to find the exact values of the control at the nodal
        % points) because is not available an analytic formula
        uMod1 = spline(tVec,uMod,tVecN0);
        uMod2 = spline(tVec,uMod,tVecN1);
        DV = ( sum(uMod1) + sum(uMod2) )*alpha;
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
    end
    
    %---------------------------------------------------------------------
    %VERIFICATION OF TRAJECTORY FROM CONTROL
    %---------------------------------------------------------------------
    if VERIFY_CONTROL    
        [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r',v');  
    end
    
end

%% 6) Calculate optimal trajectory with indirect multiple shooting

if INDIRECT_MULTIPLE_SHOOTING
    
    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  INDIRECT MULTIPLE SHOOTING %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % % SELECTION OF LAMBDA
    %--------------------------------------------------------------------------
    lambdaAd = [lambdarSS; lambdavSS];
    %lambdaAd = lambdaAd;
    
    % Define number of subintervals in which the integration is performed
    M = 80;
    
    % Set length of time vector to give to the integrator (for ease in
    % subsequent saving of r,v,lambda vectors)
    nT = 40;
    
    if SCALED_EQ   
        %--------------------------------------------------------------------------
        % VARIABLES SCALING
        %--------------------------------------------------------------------------
        % Variables scaling from AU and days 
        timeVecS = scalingFunction(timeVec,1/(tf));
        t0S = scalingFunction(t0,1/(tf));
        tfS = scalingFunction(tf,1/(tf));
        r0S = scalingFunction(r0,1/norm(r0));
        rfS = scalingFunction(rf,1/norm(r0));
        v0S = scalingFunction(v0,1/norm(v0));
        vfS = scalingFunction(vf,1/norm(v0));
        RS = scalingFunction(R,1/norm(r0));
        VcartS = scalingFunction(Vcart,1/norm(v0));
        lambdaAdS(1:3,:) = scalingFunction(lambdaAd(1:3,:),1/norm(lambdaAd(1:3,1)));
        lambdaAdS(4:6,:) = scalingFunction(lambdaAd(4:6,:),1/norm(lambdaAd(4:6,1)));
        mu = scalingFunction(mu,1/(norm(v0)^2*norm(r0)));
        
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        [r,v,lambda,tVec] = OptimalTrajectoryIndirectFSOLVEMultipleScaled(timeVecS,t0S,tfS,r0S,rfS,v0S,vfS,RS,VcartS,lambdaAdS,M,nT,mu);
        
        %--------------------------------------------------------------------------
        % VARIABLES RE-SCALING
        %--------------------------------------------------------------------------
        % Variables re-scaling into AU and days   
        tVec = scalingFunction(tVec,tf);
        r = scalingFunction(r,norm(r0));
        v = scalingFunction(v,norm(v0));
        lambda(:,1:3) = scalingFunction(lambda(:,1:3),norm(lambdaAd(1:3,1)));
        lambda(:,4:6) = scalingFunction(lambda(:,4:6),norm(lambdaAd(4:6,1)));
  
    else   
        %--------------------------------------------------------------------------
        % CHOOSE IF TO USE fmincon (Matlab) solver or OPTI
        %--------------------------------------------------------------------------  
        if OPTI
            %--------------------------------------------------------------------------
            % FUNCTION CALL
            %--------------------------------------------------------------------------
            % Call the function for indirect multiple shooting solution of optimal control
            % problem with OPTI
            [r,v,lambda,tVec] = OptimalTrajectoryIndirectOPTIMultiple(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,M,nT,mu);    
        else
            %--------------------------------------------------------------------------
            % FUNCTION CALL
            %--------------------------------------------------------------------------
            % Call the function for indirect multiple shooting solution of optimal control
            % problem with fsolve (Matlab)
            [r,v,lambda,tVec] = OptimalTrajectoryIndirectFSOLVEMultiple(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,M,nT,mu);
            
        end
        
    end
    
    %--------------------------------------------------------------------------
    % RESULTS PLOTTING
    %--------------------------------------------------------------------------
    
    % First obtain control from the optimal control condition
    Uopt = (-lambda(:,4:6))';
    
    % Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r,Uopt);
    
    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        % Transpose for compatibility
        r = r';
        v = v';
        save('traj_rvut.mat','r','v','Uopt','tVec');
        % Re-Transpose for compatibility
        r = r';
        v = v';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to obtain Uopt from the dynamics instead than from the optimal
    % control condition
    %      vDot = [];
    %      for i=1:3
    %          vDot(:,i) = ( v(2:end,i) - v(1:end-1,i) )./( tVec(2:end) - tVec(1:end-1))';
    %      end
    %      vDot(end+1,:) = ( vDot(end,:) + vDot(end-1,:) )/2;
    %
    %      Uopt = vDot + mu*r./norm(r.^3);
    %
    %      Uopt = Uopt';
    %      plotTrajectoryCartesianState(DepBody,ArrBody,r,Uopt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot comparison of r,v,u,lambda between inverse method and optimal
    % control
    comparisonLRVU(timeVec,tVec,R,r',Vcart,v',Ucart,Uopt,lambdaAd(1:3,:),lambda');
    
    %--------------------------------------------------------------------------
    % DV CALCULATION
    %--------------------------------------------------------------------------
    if CS_INT
        % Use Cavalieri-Simpson composite quadrature rule (accurate only for
        % odd number of mesh points)
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        tStep = tVec(3) - tVec(1);
        DV = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
            + 4*sum(uMod(2:2:end-1)) + uMod(end));
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
        
    else  
        % Use Gauss 2 points composite quadrature rule (accurate for odd and
        % even number of mesh points).
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        alpha = (tVec(2)-tVec(1))/2;
        tVecN0 = ( tVec(1:end-1) + tVec(2:end) )/2 + alpha/(sqrt(3));
        tVecN1 = ( tVec(1:end-1) + tVec(2:end) )/2 - alpha/(sqrt(3));
        
        % Find uMod1 and uMod2 with an interpolation (is not possible, or too
        % complicated, to find the exact values of the control at the nodal
        % points) because is not available an analytic formula
        uMod1 = spline(tVec,uMod,tVecN0);
        uMod2 = spline(tVec,uMod,tVecN1);
        DV = ( sum(uMod1) + sum(uMod2) )*alpha;
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
    end
    
    %---------------------------------------------------------------------
    %VERIFICATION OF TRAJECTORY FROM CONTROL
    %---------------------------------------------------------------------
    if VERIFY_CONTROL 
        [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r',v');
    end
    
end

%% 7) Calculate optimal trajectory with indirect multiple shooting with the 
% formulation that is linear in the control

if INDIRECT_MULTIPLE_SHOOTING_LIN_U
    
    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  INDIRECT MULTIPLE SHOOTING (H LINEAR IN U) %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % % SELECTION OF LAMBDA
    %--------------------------------------------------------------------------
    lambdaAd = [lambdarSS; lambdavSS];
    %lambdaAd = lambdaAd;
    
    % Define number of subintervals in which the integration is performed
    M = 50;
    
    % Set length of time vector to give to the integrator (for ease in
    % subsequent saving of r,v,lambda vectors)
    nT = 40;
    
    if SCALED_EQ
        %--------------------------------------------------------------------------
        % VARIABLES SCALING
        %--------------------------------------------------------------------------
        % Variables scaling from AU and days
        timeVecS = scalingFunction(timeVec,1/(tf));
        t0S = scalingFunction(t0,1/(tf));
        tfS = scalingFunction(tf,1/(tf));
        r0S = scalingFunction(r0,1/norm(r0));
        rfS = scalingFunction(rf,1/norm(r0));
        v0S = scalingFunction(v0,1/norm(v0));
        vfS = scalingFunction(vf,1/norm(v0));
        RS = scalingFunction(R,1/norm(r0));
        VcartS = scalingFunction(Vcart,1/norm(v0));
        lambdaAdS(1:3,:) = scalingFunction(lambdaAd(1:3,:),1/norm(lambdaAd(1:3,1)));
        lambdaAdS(4:6,:) = scalingFunction(lambdaAd(4:6,:),1/norm(lambdaAd(4:6,1)));
        uMaxS = scalingFunction(uMax,1/norm(lambdaAd(4:6,1)));
        %mu = scalingFunction(mu,1/(norm(v0)^2*norm(r0)));
             
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function that computes the optimal trajectory
        % Results in AU,days
        [r,v,lambda,uVec,tVec] = OptimalTrajectoryIndirectFSOLVEMultipleLinUScaled(timeVecS,t0S,tfS,r0S,rfS,v0S,vfS,RS,VcartS,lambdaAdS,M,nT,mu,uMaxS);
        
        %--------------------------------------------------------------------------
        % VARIABLES RE-SCALING
        %--------------------------------------------------------------------------
        % Variables re-scaling into AU and days  
        tVec = scalingFunction(tVec,tf);
        r = scalingFunction(r,norm(r0));
        v = scalingFunction(v,norm(v0));
        lambda(:,1:3) = scalingFunction(lambda(:,1:3),norm(lambdaAd(1:3,1)));
        lambda(:,4:6) = scalingFunction(lambda(:,4:6),norm(lambdaAd(4:6,1)));
        uVec = scalingFunction(uVec,norm(lambdaAd(4:6,1)));
             
    else 
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function that computes the optimal trajectory
        % Results in AU,days
        [r,v,lambda,uVec,tVec] = OptimalTrajectoryIndirectFSOLVEMultipleLinU(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,M,nT,mu,uMax);
        
    end
    
    %--------------------------------------------------------------------------
    % RESULTS PLOTTING
    %--------------------------------------------------------------------------
    
    % First obtain control from the control vector in output
    Uopt = uVec';
    
    % Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r,Uopt);
    
    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        % Transpose for compatibility
        r = r';
        v = v';
        save('traj_rvut.mat','r','v','Uopt','tVec');
        % Re-Transpose for compatibility
        r = r';
        v = v';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to obtain Uopt from the dynamics instead than from the optimal
    % control condition
    %      vDot = [];
    %      for i=1:3
    %          vDot(:,i) = ( v(2:end,i) - v(1:end-1,i) )./( tVec(2:end) - tVec(1:end-1))';
    %      end
    %      vDot(end+1,:) = ( vDot(end,:) + vDot(end-1,:) )/2;
    %
    %      Uopt = vDot + mu*r./norm(r.^3);
    %
    %      Uopt = Uopt';
    %      plotTrajectoryCartesianState(DepBody,ArrBody,r,Uopt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot comparison of r,v,u,lambda between inverse method and optimal
    % control
    comparisonLRVU(timeVec,tVec,R,r',Vcart,v',Ucart,Uopt,lambdaAd(1:3,:),lambda');
    
    %--------------------------------------------------------------------------
    % DV CALCULATION
    %--------------------------------------------------------------------------
    if CS_INT
        % Use Cavalieri-Simpson composite quadrature rule (accurate only for
        % odd number of mesh points)
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        tStep = tVec(3) - tVec(1);
        DV = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
            + 4*sum(uMod(2:2:end-1)) + uMod(end));
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
        
    else    
        % Use Gauss 2 points composite quadrature rule (accurate for odd and
        % even number of mesh points).
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        alpha = (tVec(2)-tVec(1))/2;
        tVecN0 = ( tVec(1:end-1) + tVec(2:end) )/2 + alpha/(sqrt(3));
        tVecN1 = ( tVec(1:end-1) + tVec(2:end) )/2 - alpha/(sqrt(3));
        
        % Find uMod1 and uMod2 with an interpolation (is not possible, or too
        % complicated, to find the exact values of the control at the nodal
        % points) because is not available an analytic formula
        uMod1 = spline(tVec,uMod,tVecN0);
        uMod2 = spline(tVec,uMod,tVecN1);
        DV = ( sum(uMod1) + sum(uMod2) )*alpha;
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
    end
    
    %---------------------------------------------------------------------
    %VERIFICATION OF TRAJECTORY FROM CONTROL
    %---------------------------------------------------------------------
    if VERIFY_CONTROL       
        [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r',v');  
    end
    
end

%% 8) Calculate optimal trajectory with solution of TPBVP with Matlab
% function bvp4c

if TPBVP
    
    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  TPBVP with Matlab bvp4c %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % SELECTION OF LAMBDA
    %--------------------------------------------------------------------------
    lambdaAd = [lambdarSS; lambdavSS];
    %lambdaAd = lambdaAd;
    
    % Define number of mesh points
    nT = 500;
    
    if SCALED_EQ
        %--------------------------------------------------------------------------
        % VARIABLES SCALING
        %--------------------------------------------------------------------------
        % Variables scaling from AU and days 
        timeVecS = scalingFunction(timeVec,1/(tf));
        t0S = scalingFunction(t0,1/(tf));
        tfS = scalingFunction(tf,1/(tf));
        r0S = scalingFunction(r0,1/norm(r0));
        rfS = scalingFunction(rf,1/norm(r0));
        v0S = scalingFunction(v0,1/norm(v0));
        vfS = scalingFunction(vf,1/norm(v0));
        RS = scalingFunction(R,1/norm(r0));
        VcartS = scalingFunction(Vcart,1/norm(v0));
        lambdaAdS(1:3,:) = scalingFunction(lambdaAd(1:3,:),1/norm(lambdaAd(1:3,1)));
        lambdaAdS(4:6,:) = scalingFunction(lambdaAd(4:6,:),1/norm(lambdaAd(4:6,1)));
        mu = scalingFunction(mu,1/(norm(v0)^2*norm(r0)));
        
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function that computes the optimal trajectory
        % Results in AU,days
        [r,v,lambda,tVec] = OptimalTrajectoryIndirectTPBVP(timeVecS,t0S,tfS,r0S,rfS,v0S,vfS,RS,VcartS,lambdaAdS,nT,mu);
        
        %--------------------------------------------------------------------------
        % VARIABLES RE-SCALING
        %--------------------------------------------------------------------------
        % Variables re-scaling into AU and days
        tVec = scalingFunction(tVec,tf);
        r = scalingFunction(r,norm(r0));
        v = scalingFunction(v,norm(v0));
        lambda(:,1:3) = scalingFunction(lambda(:,1:3),norm(lambdaAd(1:3,1)));
        lambda(:,4:6) = scalingFunction(lambda(:,4:6),norm(lambdaAd(4:6,1)));
        
    else
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function that computes the optimal trajectory
        % Results in AU,days
        [r,v,lambda,tVec] = OptimalTrajectoryIndirectTPBVP(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,nT,mu);
             
    end
    
    %--------------------------------------------------------------------------
    % RESULTS PLOTTING
    %--------------------------------------------------------------------------
    
    %First obtain control from the optimal control condition
    Uopt = (-lambda(:,4:6))';
    
    %Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r,Uopt);
    
    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        % Transpose for compatibility
        r = r';
        v = v';
        save('traj_rvut.mat','r','v','Uopt','tVec');
        % Re-Transpose for compatibility
        r = r';
        v = v';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to obtain Uopt from the dynamics instead than from the optimal
    % control condition
    %      vDot = [];
    %      for i=1:3
    %          vDot(:,i) = ( v(2:end,i) - v(1:end-1,i) )./( tVec(2:end) - tVec(1:end-1))';
    %      end
    %      vDot(end+1,:) = ( vDot(end,:) + vDot(end-1,:) )/2;
    %
    %      Uopt = vDot + mu*r./norm(r.^3);
    %
    %      Uopt = Uopt';
    %      plotTrajectoryCartesianState(DepBody,ArrBody,r,Uopt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot comparison of r,v,u,lambda between inverse method and optimal
    % control
    comparisonLRVU(timeVec,tVec,R,r',Vcart,v',Ucart,Uopt,lambdaAd(1:3,:),lambda');
    
    %--------------------------------------------------------------------------
    % DV CALCULATION
    %--------------------------------------------------------------------------
    if CS_INT
        % Use Cavalieri-Simpson composite quadrature rule (accurate only for
        % odd number of mesh points)
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        tStep = tVec(3) - tVec(1);
        DV = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
            + 4*sum(uMod(2:2:end-1)) + uMod(end));
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
        
    else   
        % Use Gauss 2 points composite quadrature rule (accurate for odd and
        % even number of mesh points).
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        alpha = (tVec(2)-tVec(1))/2;
        tVecN0 = ( tVec(1:end-1) + tVec(2:end) )/2 + alpha/(sqrt(3));
        tVecN1 = ( tVec(1:end-1) + tVec(2:end) )/2 - alpha/(sqrt(3));
        
        % Find uMod1 and uMod2 with an interpolation (is not possible, or too
        % complicated, to find the exact values of the control at the nodal
        % points) because is not available an analytic formula
        uMod1 = spline(tVec,uMod,tVecN0);
        uMod2 = spline(tVec,uMod,tVecN1);
        DV = ( sum(uMod1) + sum(uMod2) )*alpha;
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
    end
    
    %---------------------------------------------------------------------
    % VERIFICATION OF TRAJECTORY FROM CONTROL
    %---------------------------------------------------------------------
    if VERIFY_CONTROL      
     [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r',v');
    end
    
end

%% 9) Calculate optimal trajectory with direct transcription

if DIRECT_TRANSCRIPTION
    
    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  DIRECT TRANSCRIPTION %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    % Define number of points in time Vector
    N = 200;
    
    %--------------------------------------------------------------------------
    % CHOOSE IF TO USE fmincon (Matlab) solver or OPTI
    %--------------------------------------------------------------------------
    
    if OPTI
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function for direct transcritpion solution of optimal control
        % problem with OPTI
        [r,v,Uopt,tVec] = OptimalTrajectoryDirectOPTI(N,timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,Ucart,mu,uMax);
        
    else
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function for direct transcritpion solution of optimal control
        % problem with fmincon
        [r,v,Uopt,tVec] = OptimalTrajectoryDirectMat(N,timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,Ucart,mu,uMax);   
    end
    
    %--------------------------------------------------------------------------
    % RESULTS PLOTTING
    %--------------------------------------------------------------------------
    % Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r',Uopt);
    
    % Plot comparison of r,v,u,lambda between inverse method and optimal
    % control
    comparisonLRVU(timeVec,tVec,R,r,Vcart,v,Ucart,Uopt);
    
    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        save('traj_rvut.mat','r','v','Uopt','tVec');
    end
    
    %--------------------------------------------------------------------------
    % DV CALCULATION
    %--------------------------------------------------------------------------
    if CS_INT
        % Use Cavalieri-Simpson composite quadrature rule (accurate only for
        % odd number of mesh points)
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        tStep = tVec(3) - tVec(1);
        DV = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
            + 4*sum(uMod(2:2:end-1)) + uMod(end));
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
        
    else   
        % Use Gauss 2 points composite quadrature rule (accurate for odd and
        % even number of mesh points).
        uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
        alpha = (tVec(2)-tVec(1))/2;
        tVecN0 = ( tVec(1:end-1) + tVec(2:end) )/2 + alpha/(sqrt(3));
        tVecN1 = ( tVec(1:end-1) + tVec(2:end) )/2 - alpha/(sqrt(3));
        
        % Find uMod1 and uMod2 with an interpolation (is not possible, or too
        % complicated, to find the exact values of the control at the nodal
        % points) because is not available an analytic formula
        uMod1 = spline(tVec,uMod,tVecN0);
        uMod2 = spline(tVec,uMod,tVecN1);
        DV = ( sum(uMod1) + sum(uMod2) )*alpha;
        fprintf('DeltaV of the trajectory: %f km/s \n',DV*AU2km/T_sid);
    end
    
    %---------------------------------------------------------------------
    % VERIFICATION OF TRAJECTORY FROM CONTROL
    %---------------------------------------------------------------------
    if VERIFY_CONTROL     
        [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r,v);
    end
    
end

%% Generation of a tangential thrust trajectory (that does not satisfy the
% final constraints)

if TANG_THRUST_TRAJ

    %---------------------------------------------------------------------
    disp('%%%%-----------------------------------------------------%%%%');
    disp('%%%%  GENERATION OF TANGENTIAL THRUST TRAJECTORY %%%%');
    disp('%%%%-----------------------------------------------------%%%%');
    %---------------------------------------------------------------------
    
    % Integrate from r0,v0
    x0 = [r0 v0];
    nStep = 200;
    tVec = linspace(t0,tf,nStep);
    options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [tOde,YOde] = ode45(@(t,x) diffEqTangU(t,x,mu,uMax), tVec, x0);
    
    % Extract position and velocity
    r = (YOde(:,1:3))';
    v = (YOde(:,4:6))';
    
    % Find the control
    for i=1:length(tVec)
        normV = norm(v(:,i));
        unitV = v(:,i)/normV;
        Uopt(:,i) = uMax*unitV;
    end
    
    % Plot trajectory obtained
    plotTrajectoryCartesianState(DepBody,ArrBody,t0,TOF,r',Uopt);
    
    %--------------------------------------------------------------------------
    % SAVING OF VARIABLES
    %--------------------------------------------------------------------------
    if SAVE_TRAJ
        save('traj_rvut.mat','r','v','Uopt','tVec');
    end
    
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end