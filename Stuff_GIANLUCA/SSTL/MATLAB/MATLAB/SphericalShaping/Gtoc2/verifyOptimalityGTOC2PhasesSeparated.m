function verifyOptimalityGTOC2PhasesSeparated(DepBody,ArrBody,RMatLegs,VcartLegs,UcartLegs,timeLegs,mu,DVLegs,depDate_vec,arrDate_vec)
%
% verifyOptimalityGTOC2PhasesSeparated
% Version for the GTOC2 asteroids randez-vous
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
% The methods used are indirect multiple shooting and direct transcription
% The different phases are solved for separately
%
% INPUT
% DepBody, ArrBody: departure and arrival bodies
% RMatLegs: multidimensional matrix containing the radius vector for each
% leg
% VcartLegs: multidimensional matrix containing the radius vector for each
% leg
% UcartLegs: multidimensional matrix containing the radius vector for each
% legtimeLegs
% timeLegs: matrix containing the time vector for each phase
% mu: gravitational parameter
% DVLegs: vector containign the DV for each leg
% depDate_vec,arrDate_vec: vector of arrival and departure dates
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
INDIRECT_MULTIPLE_SHOOTING = 0;
DIRECT_TRANSCRIPTION = 0;

% Select the Optimal Control Method
choice = menu('Select Optimal Control Method','INDIRECT_MULTIPLE_SHOOTING','DIRECT_TRANSCRIPTION','NO METHOD');
switch choice
    case 1
        INDIRECT_MULTIPLE_SHOOTING = 1;
    case 2
        DIRECT_TRANSCRIPTION = 1;
    case 3
end

%--------------------------------------------------------------------------
% FLAG SELECTION
%--------------------------------------------------------------------------
% Select (1) or not select (0) if to:
% 1) Integrate forward or backward the diff eq
% 2) Integrate only the adjoints diff eq or the complete diff eq
% 3) Plot the adjoints and Hamiltonian
% 4) Verify the control obtained by the selected direct or indirect method
% 5) Quadrature method for DV: 1 for CS composite, 0 for Gauss composite
% 6) Use the alternative OPTI Toolbox for solving (where implemented)
FORWARD_INT = 0;
ONLY_ADJOINTS = 1;
PLOT_ADJOINTS_HAM = 0;
VERIFY_CONTROL = 1;
CS_INT = 0;
OPTI = 0;

%% 2) Calculate the history of the Langrange multipliers from the control
% obtained (multipliers correspondent to the shaped solution)

%--------------------------------------------------------------------------
% Cycle on the different legs
%--------------------------------------------------------------------------
for Nleg=1:4
    
    % Write ausiliary variables for R,Vcart,Ucart
    R_aux(:,:) = RMatLegs(Nleg,:,:);
    V_aux(:,:) = VcartLegs(Nleg,:,:);
    U_aux(:,:) = UcartLegs(Nleg,:,:);
    
    lambdavSS(Nleg,:,:) = -UcartLegs(Nleg,:,:);
    lambdavdotSS = [];
    
    % Cycle to approximate the derivative of lambdav with finite
    % differences
    for i=1:3
        DeltaLambdaV(1,:) = lambdavSS(Nleg,i,2:end) - lambdavSS(Nleg,i,1:end-1);
        lambdavdotSS(Nleg,i,:) = ( DeltaLambdaV )./( timeLegs(Nleg,2:end) - timeLegs(Nleg,1:end-1));
    end
    
    % Add last element (not available from finite differencing
    % Use interpolation to extrapolate (doesnt work for now)
    %lambdavdotSS(:,end+1) = interp1(timeVec(1:end-1),lambdavdotSS(:,:),timeVec(end),'nearest','extrap');
    % Approximate as a mean of the previous two derivatives (very rough)
    lambdavdotSS(Nleg,:,end+1) = ( lambdavdotSS(Nleg,:,end) + lambdavdotSS(Nleg,:,end-1) )/2;
    
    % Obtain lambda r
    lambdarSS(Nleg,:,:) = -lambdavdotSS(Nleg,:,:);
    
    % Auxiliary lambdar,lambdav
    lambdar_aux(:,:) = lambdarSS(Nleg,:,:);
    lambdav_aux(:,:) = lambdavSS(Nleg,:,:);
    
%% 3) Calculate the history of the Lagrange multipliers from the adjoint differential equations
    
    % Choose starting conditions depending upon integration direction
    if (FORWARD_INT == 1)
        
        % Initial conditions for integration forward
        lambdar0(:,1) = lambdarSS(Nleg,:,1);
        lambdav0(:,1) = lambdavSS(Nleg,:,1);
        
        r0(:,1) = RMatLegs(Nleg,:,1);
        v0(:,1) = VcartLegs(Nleg,:,1);
        
        % Create the Flip vector (even if in this case is not flipped)
        timeVecFlip(1,:) = timeLegs(Nleg,:);
        
    else
        
        % Initial conditions for integration backward
        lambdar0(:,1) = lambdarSS(Nleg,:,end);
        lambdav0(:,1) = lambdavSS(Nleg,:,end);
        
        r0(:,1) = RMatLegs(Nleg,:,end);
        v0(:,1) = VcartLegs(Nleg,:,end);
        
        % Flip the time vector for backward integration
        timeVecFlip(1,:) = fliplr(timeLegs(Nleg,:));
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
        [tOde,lambdaAdOde] = ode45(@(t,lambda) adjointEqOde(t,lambda,R_aux,timeLegs(Nleg,:),mu), timeVecFlip, x0);
        
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
            plot(timeLegs(Nleg,:),R_aux(1,:),'b','linewidth',2);
            hold on;
            plot(t,Y(1,:),'g','linewidth',2);
            ylabel('R x');
            legend('Inverse Method','Int From Adjoints','location','northeast');
            subplot(3,1,2);
            plot(timeLegs(Nleg,:),R_aux(2,:),'b','linewidth',2);
            hold on;
            plot(t,Y(2,:),'g','linewidth',2);
            ylabel('R y');
            subplot(3,1,3);
            plot(timeLegs(Nleg,:),R_aux(3,:),'b','linewidth',2);
            hold on;
            plot(t,Y(3,:),'g','linewidth',2);
            ylabel('R z');
            tit = strcat('Radius vector comparison Inverse Method-Integration From Adjoints for Leg N',int2str(Nleg));
            suptitle(tit);
            xlabel('Time');
            
            %Radius vector
            figure;
            subplot(3,1,1);
            plot(timeLegs(Nleg,:),V_aux(1,:),'b','linewidth',2);
            hold on;
            plot(t,Y(4,:),'g','linewidth',2);
            ylabel('V x');
            legend('Inverse Method','Int From Adjoints','location','northeast');
            subplot(3,1,2);
            plot(timeLegs(Nleg,:),V_aux(1,:),'b','linewidth',2);
            hold on;
            plot(t,Y(5,:),'g','linewidth',2);
            ylabel('V y');
            subplot(3,1,3);
            plot(timeLegs(Nleg,:),V_aux(1,:),'b','linewidth',2);
            hold on;
            plot(t,Y(6,:),'g','linewidth',2);
            ylabel('V z');
            tit = strcat('Velocity vector comparison Inverse Method-Integration From Adjoints for Leg N',int2str(Nleg));
            suptitle(tit);
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
        plot(timeLegs(Nleg,:),lambdar_aux(1,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(1,:),'g','linewidth',2);
        ylabel('\lambda rx');
        legend('Opt control condition','Adjoint diff equations','location','northeast');
        subplot(3,1,2);
        plot(timeLegs(Nleg,:),lambdar_aux(2,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(2,:),'g','linewidth',2);
        ylabel('\lambda ry');
        subplot(3,1,3);
        plot(timeLegs(Nleg,:),lambdar_aux(3,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(3,:),'g','linewidth',2);
        ylabel('\lambda rz');
        tit = strcat('r adjoints comparison Inverse Method-Integration From Adjoints for Leg N',int2str(Nleg));
        suptitle(tit);
        xlabel('Time');
        
        % V adjoints
        figure;
        subplot(3,1,1);
        plot(timeLegs(Nleg,:),lambdav_aux(1,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(4,:),'g','linewidth',2);
        ylabel('\lambda vx');
        legend('Opt control condition','Adjoint diff equations','location','northeast');
        subplot(3,1,2);
        plot(timeLegs(Nleg,:),lambdav_aux(2,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(5,:),'g','linewidth',2);
        ylabel('\lambda vy');
        subplot(3,1,3);
        plot(timeLegs(Nleg,:),lambdav_aux(3,:),'b','linewidth',2);
        hold on;
        plot(t,lambdaAd(6,:),'g','linewidth',2);
        ylabel('\lambda vz');
        tit = strcat('v adjoints comparison Inverse Method-Integration From Adjoints for Leg N',int2str(Nleg));
        suptitle(tit);
        xlabel('Time');
        
        %--------------------------------------------------------------------------
        % HAMILTONIAN CALCULATION
        %--------------------------------------------------------------------------
        %Calculate the Hamiltonian at each moment of time
        % with the adjoints coming from the optimal control condition
        HSS = [];
        HAd = [];
        for i=1:length(timeLegs(Nleg,:))
            HSS(i) = 0.5*norm(U_aux(:,i))^2 + dot( lambdar_aux(:,i),V_aux(:,i) ) + dot( lambdav_aux(:,i), (  -mu*R_aux(:,i)/norm(R_aux(:,i))^3  + U_aux(:,i)) );
            HAd(i) = 0.5*norm(U_aux(:,i))^2 + dot( lambdaAd(1:3,i),V_aux(:,i))  + dot( lambdaAd(4:6,i), (  -mu*R_aux(:,i)/norm(R_aux(:,i))^3  + U_aux(:,i)) );
        end
        
        % Hamiltonian plotting
        figure;
        plot(timeLegs(Nleg,:),HSS,'b','linewidth',2);
        hold on;
        plot(t,HAd,'g','linewidth',2);
        plot(timeLegs(Nleg,:),ones(1,length(timeLegs(Nleg,:)))*HSS(1),'r','linewidth',2);
        xlabel('Time');
        ylabel('Hamiltonian function');
        tit = strcat('Hamiltonian Leg N',int2str(Nleg));
        suptitle(tit);
        legend('H inv','H adj','Initial value');
        
    end
    
end

%% Set BC and parameters for Optimal Trajectory solution

for Nleg=1:4
    %--------------------------------------------------------------------------
    % BOUNDARY CONDITIONS AND PARAMETERS
    %--------------------------------------------------------------------------
    % Write ausiliary variables for R,Vcart,Ucart,lambdar,lambdav
    R_aux(:,:) = RMatLegs(Nleg,:,:);
    V_aux(:,:) = VcartLegs(Nleg,:,:);
    U_aux(:,:) = UcartLegs(Nleg,:,:);
    lambdar_aux(:,:) = lambdarSS(Nleg,:,:);
    lambdav_aux(:,:) = lambdavSS(Nleg,:,:);
    
    % Write boundary conditions
    % To be more correct, the initial and final positions and velocitites
    % should be calculated from the Departure and Arrival bodies, but since
    % the inverse solution has been verified, it is known that such boundaries
    % are exact (within tolerances)
    t0 = timeLegs(Nleg,1);
    tf = timeLegs(Nleg,end);
    TOF = tf - t0;
    r0 = R_aux(:,1)';
    v0 = V_aux(:,1)';
    rf = R_aux(:,end)';
    vf = V_aux(:,end)';
    uMax = DVLegs(Nleg)/TOF;
    
    %% 6) Calculate optimal trajectory with indirect multiple shooting
    
    if INDIRECT_MULTIPLE_SHOOTING
        
        %---------------------------------------------------------------------
        disp('%%%%-----------------------------------------------------%%%%');
        disp('%%%%  INDIRECT MULTIPLE SHOOTING %%%%');
        disp('%%%%-----------------------------------------------------%%%%');
        %---------------------------------------------------------------------
        fprintf('Search for Leg N %g \n',Nleg);
        disp('%%%%-----------------------------------------------------%%%%');
        %---------------------------------------------------------------------
        
        %--------------------------------------------------------------------------
        % % SELECTION OF LAMBDA
        %--------------------------------------------------------------------------
        lambdaAd = [lambdar_aux; lambdav_aux];
        %lambdaAd = lambdaAd;
        
        % Define number of subintervals in which the integration is performed
        M = 60;
        
        % Set length of time vector to give to the integrator (for ease in
        % subsequent saving of r,v,lambda vectors)
        nT = 40;
        
        %--------------------------------------------------------------------------
        % FUNCTION CALL
        %--------------------------------------------------------------------------
        % Call the function for indirect multiple shooting solution of optimal control
        % problem with fsolve (Matlab)
        [r,v,lambda,tVec,exitFlag] = OptimalTrajectoryIndirectFSOLVEMultiple(timeLegs(Nleg,:),t0,tf,r0,rf,v0,vf,R_aux,V_aux,lambdaAd,M,nT,mu);
              
        %--------------------------------------------------------------------------
        % DV CALCULATION
        %--------------------------------------------------------------------------
        
        % First obtain control from the optimal control condition
        Uopt = (-lambda(:,4:6))';
        
        if CS_INT
            % Use Cavalieri-Simpson composite quadrature rule (accurate only for
            % odd number of mesh points)
            uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
            tStep = tVec(3) - tVec(1);
            DVNewLegs(Nleg) = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
                + 4*sum(uMod(2:2:end-1)) + uMod(end));
            
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
            DVNewLegs(Nleg) = ( sum(uMod1) + sum(uMod2) )*alpha;

        end
        
        %--------------------------------------------------------------------------
        % RESULTS PLOTTING AND VERIFICATION OF TRAJECTORY FROM CONTROL
        %--------------------------------------------------------------------------
        
        % Plot comparison of r,v,u,lambda between inverse method and optimal
        % control
        comparisonLRVU(timeLegs(Nleg,:),tVec,R_aux,r',V_aux,v',U_aux,Uopt,lambdaAd(1:3,:),lambda');
        
        if VERIFY_CONTROL
            TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r',v');
        end
        
        %--------------------------------------------------------------------------
        % SAVE RESULTS IN MULTIDIMENSIONAL MATRIX
        %--------------------------------------------------------------------------
        R_multi(Nleg,:,:) = r';
        V_multi(Nleg,:,:) = v';
        U_multi(Nleg,:,:) = Uopt;
        time_multi(Nleg,:) = tVec;
        exitFlag_multi(Nleg) = exitFlag;
    end
    
    
    %% 9) Calculate optimal trajectory with direct transcription
    
    if DIRECT_TRANSCRIPTION
        
        %---------------------------------------------------------------------
        disp('%%%%-----------------------------------------------------%%%%');
        disp('%%%%  DIRECT TRANSCRIPTION %%%%');
        disp('%%%%-----------------------------------------------------%%%%');
        %---------------------------------------------------------------------
        %---------------------------------------------------------------------
        fprintf('Search for Leg N %g \n',Nleg);
        disp('%%%%-----------------------------------------------------%%%%');
        %---------------------------------------------------------------------
        
        % Define number of points in time Vector
        N = 150;
        
        %--------------------------------------------------------------------------
        % CHOOSE IF TO USE fmincon (Matlab) solver or OPTI
        %--------------------------------------------------------------------------  
        if OPTI
            %--------------------------------------------------------------------------
            % FUNCTION CALL
            %--------------------------------------------------------------------------
            % Call the function for direct transcritpion solution of optimal control
            % problem with OPTI
            [r,v,Uopt,tVec,exitFlag] = OptimalTrajectoryDirectOPTI(N,timeLegs(Nleg,:),t0,tf,r0,rf,v0,vf,R_aux,V_aux,U_aux,mu,uMax);
            
        else
            %--------------------------------------------------------------------------
            % FUNCTION CALL
            %--------------------------------------------------------------------------
            % Call the function for direct transcritpion solution of optimal control
            % problem with fmincon
            [r,v,Uopt,tVec,exitFlag] = OptimalTrajectoryDirectMat(N,timeLegs(Nleg,:),t0,tf,r0,rf,v0,vf,R_aux,V_aux,U_aux,mu,uMax);
        end
             
        %--------------------------------------------------------------------------
        % DV CALCULATION
        %--------------------------------------------------------------------------
        if CS_INT
            % Use Cavalieri-Simpson composite quadrature rule (accurate only for
            % odd number of mesh points)
            uMod = sqrt(Uopt(1,:).^2 + Uopt(2,:).^2 + Uopt(3,:).^2);
            tStep = tVec(3) - tVec(1);
            DVNewLegs(Nleg) = (tStep/6)*(uMod(1) + 2*sum( uMod(3:2:end-1) ) + ...
                + 4*sum(uMod(2:2:end-1)) + uMod(end));
            
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
            DVNewLegs(Nleg) = ( sum(uMod1) + sum(uMod2) )*alpha;

        end
        
        %---------------------------------------------------------------------
        % RESULTS PLOTTING AND VERIFICATION OF TRAJECTORY FROM CONTROL
        %---------------------------------------------------------------------
        
        % Plot comparison of r,v,u,lambda between inverse method and optimal
        % control
        comparisonLRVU(timeLegs(Nleg,:),tVec,R_aux,r,V_aux,v,U_aux,Uopt);
        
        if VERIFY_CONTROL
            TrajectoryFromControl(tVec,r0,v0,Uopt,mu,r,v);
        end
        
        %--------------------------------------------------------------------------
        % SAVE RESULTS IN MULTIDIMENSIONAL MATRIX
        %--------------------------------------------------------------------------
        R_multi(Nleg,:,:) = r;
        V_multi(Nleg,:,:) = v;
        U_multi(Nleg,:,:) = Uopt;
        time_multi(Nleg,:) = tVec;
        exitFlag_multi(Nleg) = exitFlag;
    end
    

end

%% Plot the final trajectory obtained

% Plot trajectory
plotTrajectoryMultBodiesCartContMultTraj(DepBody,ArrBody,depDate_vec,arrDate_vec,R_multi,U_multi);

% Show information on convergence for each leg
disp('%%%%-----------------------------------------------------%%%%');
for i=1:4
 fprintf('Optimization on Leg N %g : Exit for %g \n',i,exitFlag_multi(i));
 fprintf('DeltaV of the Leg N %g is %f km/s \n',i,DVNewLegs(i)*AU2km/T_sid);
 disp('%%%%%%%%%%%%%%%%%%%%%%%%%');
end
disp('%%%%-----------------------------------------------------%%%%');
fprintf('Total DeltaV is %f km/s \n',sum(DVNewLegs)*AU2km/T_sid);
disp('%%%%-----------------------------------------------------%%%%');

end