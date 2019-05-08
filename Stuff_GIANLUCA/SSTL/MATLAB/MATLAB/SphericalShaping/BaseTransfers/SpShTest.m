%Spherical Shaping Test
%28/09/15
%-------------------------------------------------------------------------%
% Apply the Spherical Shaping method to the test cases of the article
% (Novak,Vasile-2011)
% First load the data, then select the target between the ones
% available, then choose other relevant parameters (time dep, time of
% flight (TOF), number of revs (nr), tolerances, angle step, note: the param
% maxIter is meaningful only if the version with manual Netwton-Raphson 
% cycle is used).
% Calculate the trajectory calling the spherical shaping function.
% At the end plot the trajectory and the graph with the DeltaV.
% In the second part of the file there is the function that verifies 
% the trajectory computed and obtain the optimal trajectory (choosing 
% between a certain amount of method).
% In the last part there is a grid search on the transfer Earth-selected 
% celestial body (with variable dep time, TOF, nr).

%NOTE: due to scaling purposes, the units used for the computation are AU
%and days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and parameters
clc;
clear all;
close all;

% Add folders containing subfunctions
%addpath(genpath('SphericalShapingFunctions'));
addpath('SphericalShapingFunctions');
addpath('TrajectoryPlotting');
addpath('OptimalityFunctions');

% Load the data
LoadData;

% Set the parameters of the spacecraft
m_i = 1000;                  % kg
Ceff = 3000*9.81*1e-3;       % Ceff = Isp*g0 [km/s]

% Set tolerance,step and max iterations
toll = 1e-3;                % at least 1e-3 (1 day of error for TOF=100d)
thetaStep = 2*pi/200;
maxIter = 600;

% Set the departure time, TOF (MJD2000 and days), and num revs
t_dep = 8000;
TOF = 900;
nr = 1;

% Set coefficient for Newton-Rapshon (if the version of the algorithm with
% Newton-Raphson is used)
a2 = 0.0;
h = 0.001;
plotIter = 0;

%-------------------------------------------------------------------------%
% Select Departure and Arrival Bodies (to use other bodies add them in the
% file LoadData)
% Define Earth as departure body
DepBody = body(1);

% Select the target body
choice = menu('Select the target body','Mars','1989ML','Tempel-1','Neptune');
switch choice
    case 1
       ArrBody = body(2); 
    case 2
       ArrBody = body(3);
    case 3
       ArrBody = body(4); 
    case 4
       ArrBody = body(5); 
end

%% Call the spherical shaping function (for a fixed t_dep and TOF)
disp('SPHERICAL SHAPING METHOD TEST CASES');
disp('-------------------------------------------------------------------');
disp(['TRIAL ON A SINGLE TRANSFER: ',DepBody.name,' - ', ArrBody.name]);
fprintf('Dep date \t \t TOF \t \t Num revs \t \t a2_0 \n');
fprintf('%f \t %f \t %g \t \t %f \n',t_dep,TOF,nr,a2);
disp('-------------------------------------------------------------------');

%[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingNoStepAcc(DepBody,ArrBody,t_dep,TOF,nr,toll,thetaStep,maxIter,a2,h,plotIter);
%[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingWithFMincon(DepBody,ArrBody,t_dep,TOF,nr,toll,thetaStep,a2);
tic
[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingWithFSolve(DepBody,ArrBody,t_dep,TOF,nr,toll,thetaStep,a2);
toc

% Check if the method is converged
if isLoopConverged==0
    fprintf('--------------------------------------------------------------- \n');
    disp('ERROR - The loop has not converged; no trajectory found with the desired parameters');
    fprintf('The number of iterations done is: %g / %g   \n',cont,maxIter); 
    fprintf('--------------------------------------------------------------- \n');
    
else
    fprintf('--------------------------------------------------------------- \n');
    disp('CONVERGENCE - The loop has converged; a trajectory has been found');
    fprintf('The number of iterations done is: %g / %g   \n',cont,maxIter); 
    
    % Calculate the final mass
    m_f = m_i*exp(-DV*AU2km/T_sid/Ceff);
    
    % Print data
    fprintf('The DeltaV required and the final mass are \n');
    fprintf('DeltaV: %f  km/s \n',DV*AU2km/T_sid); 
    fprintf('Final mass: %f  kg \n',m_f);
    fprintf('--------------------------------------------------------------- \n');
    
    % Plot the trajectory for Departure and Arrival bodies
    plotTrajectory(DepBody,ArrBody,t_dep,TOF,thetaVec,R,Phi,Ucart);
    %plotTrajectoryMovie(DepBody,ArrBody,thetaVec,R,Phi,mu_S,t_dep,TOF,Ucart);
    
    % Plot the control as a function of angle theta (not time!!)
    figure;
    subplot(3,1,1);
    plot(thetaVec*r2d,u_mat(1,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Ut [km/s^2]');
    subplot(3,1,2);
    plot(thetaVec*r2d,u_mat(2,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Un [km/s^2]');
    subplot(3,1,3);
    plot(thetaVec*r2d,u_mat(3,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uh [km/s^2]');
    xlabel('Theta [deg]');
    suptitle('Control acceleration');
    set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
    
%    Plot the control as a function of angle theta (not time!!) in
%   cartesian coordinates
    figure;
    subplot(3,1,1);
    plot(thetaVec*r2d,Ucart(1,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Ux [km/s^2]');
    subplot(3,1,2);
    plot(thetaVec*r2d,Ucart(2,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uy [km/s^2]');
    subplot(3,1,3);
    plot(thetaVec*r2d,Ucart(3,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uz [km/s^2]');
    xlabel('Theta [deg]');
    suptitle('Control acceleration in cart coordinates');
    set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');

    % Plot the modulus of the control acceleration
    figure;
    plot(thetaVec*r2d,sqrt(Ucart(1,:).^2 + Ucart(2,:).^2 + Ucart(3,:).^2)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('U [km/s^2]');
    xlabel('Theta [deg]');
    suptitle('Modulus of control acceleration from Spherical Shaping');

keyboard
% (IF OPTIMALITY VERIFICATION AND COMPUTATION
%      IS NOT WANTED, STOP CODE WITH KEYBOARD) 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verify optimality of the obtained shaped solution respect to optimal
    % control equations and compute optimal trajectory

    % First obtain the time history vector
    [timeVec] = timeEvolution(thetaVec,TPrime,t_dep,TOF,toll);
    
    % Obtain trajectory as cartesian state (x,y,z)
    Xt = R.*cos(Phi).*cos(thetaVec);
    Yt = R.*cos(Phi).*sin(thetaVec);
    Zt = R.*sin(Phi);
    
    Rmat = [];
    Rmat(1,:) = Xt;
    Rmat(2,:) = Yt;
    Rmat(3,:) = Zt;   
    
    verifyOptimality(DepBody,ArrBody,Rmat,Vcart,Ucart,timeVec,mu_S,DV);
    
    %---------------------------------------------------------------------
    % Verify the spherical shaping solution
     
%      tVec = linspace(timeVec(1),timeVec(end),length(timeVec));
%      Ucart1(1,:) = spline(timeVec,Ucart(1,:),tVec);
%      Ucart1(2,:) = spline(timeVec,Ucart(2,:),tVec);
%      Ucart1(3,:) = spline(timeVec,Ucart(3,:),tVec);
%      r0 = Rmat(:,1)';
%      v0 = Vcart(:,1)';
%      [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Ucart1,mu_S);
%      
%      % Radius vector
%      figure;
%      subplot(3,1,1);
%      plot(timeVec,Rmat(1,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,R_ver(1,:),'g','linewidth',2);
%      ylabel('R x');
%      legend('Optimal Control','IntWControl','location','northeast');
%      subplot(3,1,2);
%      plot(timeVec,Rmat(2,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,R_ver(2,:),'g','linewidth',2);
%      ylabel('R y');
%      subplot(3,1,3);
%      plot(timeVec,Rmat(3,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,R_ver(3,:),'g','linewidth',2);
%      ylabel('R z');
%      suptitle('Radius vector comparison Sp Shaping-Integration with Control');
%      xlabel('Time');
%           
%      % Velocity vector
%      figure;
%      subplot(3,1,1);
%      plot(timeVec,Vcart(1,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,V_ver(1,:),'g','linewidth',2);
%      ylabel('V x');
%      legend('Optimal Control','IntWControl','location','northeast');
%      subplot(3,1,2);
%      plot(timeVec,Vcart(2,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,V_ver(2,:),'g','linewidth',2);
%      ylabel('V y');
%      subplot(3,1,3);
%      plot(timeVec,Vcart(3,:),'b','linewidth',2);
%      hold on;
%      plot(tVec,V_ver(3,:),'g','linewidth',2);
%      ylabel('V z');
%      suptitle('Velocity vector comparison Sp Shaping-Integration with Control');
%      xlabel('Time');
    %---------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end

keyboard
% (IF GRID SEARCH IS NOT WANTED, STOP CODE WITH KEYBOARD)

%% Cycle on t_dep,TOF,nr to obtain the plot of DeltaV respect to departure data and time of flight
disp(['GRID SEARCH ON THE TRANSFER: ',DepBody.name,' - ', ArrBody.name]);
fprintf('--------------------------------------------------------------- \n');
%fprintf('Percentage completed \t');

% Set vectors of DepTime and TOF
t_dep_vect = [7500:40:10000];     %MJD2000 [7305:50:10226]
TOF_vect = [11000:200:30000];        %days    [500:40:2000];
Ndep = length(t_dep_vect);
Ntof = length(TOF_vect);

% Set vectors of number of revolutions
Nr_vect = [0:1];

% Set tolerance,step and max iterations
toll = 1e-3;                % at least 1e-3 (1 day of error for TOF=1000d)
thetaStep = 2*pi/100;
maxIter = 500;

% Set coefficient for Newton-Rapshon (if the version of the algorithm with
% Newton-Raphson is used)
a2 = 0.0;
h = 0.001;
plotIter = 0;

% Initialize other variables
convPerc = 0;
DV_matrix = [];

% Start the cycle
tic
for i=1:Ntof
    for j=1:Ndep
        DV_min = 1e6;
        convergenceForOneNr = 0;
        for k=1:length(Nr_vect)
            %[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingNoStepAcc(DepBody,ArrBody,t_dep_vect(j),TOF_vect(i),Nr_vect(k),toll,thetaStep,maxIter,a2,h,plotIter);
            %[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingWithFMincon(DepBody,ArrBody,t_dep_vect(j),TOF_vect(i),Nr_vect(k),toll,thetaStep,a2);
            [thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingWithFSolve(DepBody,ArrBody,t_dep_vect(j),TOF_vect(i),Nr_vect(k),toll,thetaStep,a2);
            if isLoopConverged==1
                if DV < DV_min
                    convergenceForOneNr = 1;
                    DV_min = DV;
                end
            end
        end
        if convergenceForOneNr==1   %If converged for at least one Nr
            convPerc = convPerc+1;
            DV_matrix(i,j) = DV_min;
        else
            % If not converged, then set NaN so that the correct minimum is
            % found and the plot doesnt show singularities (NaN is plotted in
            % white)
            DV_matrix(i,j) = nan;
        end
    end
    % Show the percentage (but only each Ndep internal loops)
    fprintf(' %f %% \n',i/Ntof*100);
end

TIMEcomp = toc;
fprintf('Results: \n');
fprintf('Total computational time: %f min \n',TIMEcomp/60);
fprintf('Comp time per trajectory: %f s \n',TIMEcomp/(Ntof*Ndep));
fprintf('Number of trajectories: %f, / %f\n', convPerc, Ntof*Ndep);
fprintf('Percentage of feasible trajectories: %f %% \n',100*convPerc/(Ntof*Ndep));
fprintf('DeltaV of best trajectory: %f km/s \n',min(min(DV_matrix))*AU2km/T_sid);
fprintf('--------------------------------------------------------------- \n');

% Plot the graphs
for i=1:Ntof
   X(i,:) = t_dep_vect; 
end
for i=1:Ndep
   Y(:,i) = TOF_vect; 
end
figure;
contourf(X,Y,DV_matrix*AU2km/T_sid);
c = colorbar;
title('DeltaV obtained [km/s]');
xlabel('Departure date [MJD2000]');
ylabel('Transfer time [days]');
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');

%-------------------------------------------------------------------------%
% Filter out DV higher than 50 km/s
for i=1:Ntof
    for j=1:Ndep
        if DV_matrix(i,j)*AU2km/T_sid > 50
           DV_matrix(i,j) = nan;
        end
    end
end

% Plot the graphs
figure;
contourf(X,Y,log(DV_matrix*AU2km/T_sid));
%c = colorbar;
title('DeltaV obtained [km/s] constrained to be under 50 km/s');
xlabel('Departure date [MJD2000]');
ylabel('Transfer time [days]');
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');