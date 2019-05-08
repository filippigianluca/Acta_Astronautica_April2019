% Test_Hill_Shaping_1: test on the shaping method for the Pseudo Equinoctal Parameters
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and parameters
clc;
clear all;
close all;

% Load the data
LoadData;

% Set the parameters of the spacecraft
m_i = 1000;                  % kg
Ceff = 3000*9.81*1e-3;       % Ceff = Isp*g0 [km/s]

% Set tolerance,step and max iterations
toll = 1e-3;                % at least 1e-3 (1 day of error for TOF=100d)
LStep = 2*pi/200;
maxIter = 100;

% Set the departure time, TOF (MJD2000 and days), and num revs
t_dep = 8000;
TOF = 900;
nr = 1;
timeStep = 1;

% Choose if analytic computation of V,A and type of shape and the lambda0
Analytic_VA = 0;
Lin_Trig = 1;
lambda0 = [0.002 0.002 0.002];

%% Bodies selection

%-------------------------------------------------------------------------
% Departure and Arrival bodies
%-------------------------------------------------------------------------
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


%% Execute the shaping

% Find departure and arrival position

%--------------------------------------------------------------------------
% Departure Body
%--------------------------------------------------------------------------
% Position and velocity of Dep Body
M_dep = DepBody.M0 + DepBody.n*(t_dep-DepBody.t0);
Theta_dep = M2theta(M_dep,DepBody.e);

% Obtain Hill Parameters from Kepler Parameters
Pseudo_Dep = Kep2Pseudo([DepBody.a DepBody.e DepBody.i DepBody.Omega DepBody.w Theta_dep]);

%--------------------------------------------------------------------------
% Arrival Body
%--------------------------------------------------------------------------
% Position and velocity of Dep Body
M_arr = ArrBody.M0 + ArrBody.n*(t_dep+TOF-ArrBody.t0);
Theta_arr = M2theta(M_arr,ArrBody.e);

% Obtain Hill Parameters from Kepler Parameters
Pseudo_Arr = Kep2Pseudo([ArrBody.a ArrBody.e ArrBody.i ArrBody.Omega ArrBody.w Theta_arr]);

%--------------------------------------------------------------------------
% Function call
%--------------------------------------------------------------------------
% Choose the shape to use (Lin_Trig or Exponential)
[LVec,timeVec,Pseudo_Matrix,R,V,ACont,DV,isLoopConverged] = Pseudo_Shaping_FMINCON(mu_S,Pseudo_Dep,Pseudo_Arr,t_dep,TOF,nr,toll,LStep,lambda0,Lin_Trig,Analytic_VA);

%% Results plotting

if isLoopConverged
    
    % Plot trajectory
    plotTrajectoryCartesianState(DepBody,ArrBody,t_dep,TOF,R,ACont');
    
    % Calculate the final mass
    m_f = m_i*exp(-DV*AU2km/T_sid/Ceff);
    
    % Print data
    fprintf('The DeltaV required and the final mass are \n');
    fprintf('DeltaV: %f  km/s \n',DV*AU2km/T_sid);
    fprintf('Final mass: %f  kg \n',m_f);
    fprintf('--------------------------------------------------------------- \n');
 
    %--------------------------------------------------------------------------
    % Verify the trajectory with the control obtained
    %--------------------------------------------------------------------------
    % Before interpolate
%     tVec = linspace(timeVec(1),timeVec(end),length(timeVec));
%     Ucart(1,:) = spline(timeVec,ACont(:,1),tVec);
%     Ucart(2,:) = spline(timeVec,ACont(:,2),tVec);
%     Ucart(3,:) = spline(timeVec,ACont(:,3),tVec);

    r0 = R(1,:);
    v0 = V(1,:);
    [R_ver,V_ver] = TrajectoryFromControl(timeVec,r0,v0,ACont',mu_S,R',V');
    
    plotTrajectoryCartesianState(DepBody,ArrBody,t_dep,TOF,R_ver',ACont');
    
    %--------------------------------------------------------------------------
    % Plot control in tangential,normal,out of plane direction
    %--------------------------------------------------------------------------
    
    % Find at each moment of time the t-n-h directions
    for i=1:length(timeVec)
       dirT(i,:) = V(i,:)/norm(V(i,:)); 
       dirR(i,:) = R(i,:)/norm(R(i,:));
       dirH(i,:) = cross(dirT(i,:),dirR(i,:));
       dirH(i,:) = dirH(i,:)/norm(dirH(i,:));
       dirN(i,:) = cross(dirH(i,:),dirT(i,:));
       dirN(i,:) = dirN(i,:)/norm(dirN(i,:));
       
       % Project Acceleration on the versors
       A_TNH(i,1) = dot(ACont(i,:),dirT(i,:));
       A_TNH(i,2) = dot(ACont(i,:),dirN(i,:));
       A_TNH(i,3) = dot(ACont(i,:),dirH(i,:));
    end
    
    figure;
    subplot(3,1,1);
    plot(timeVec,A_TNH(:,1)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Ut [km/s^2]');
    subplot(3,1,2);
    plot(timeVec,A_TNH(:,2)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Un [km/s^2]');
    subplot(3,1,3);
    plot(timeVec,A_TNH(:,3)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uh [km/s^2]');
    xlabel('Time [days]');
    suptitle('Control acceleration in t-n-h coordinates');
    set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
    
    %--------------------------------------------------------------------------
    % Plot control in cartesian coordinates
    %--------------------------------------------------------------------------
    figure;
    subplot(3,1,1);
    plot(timeVec,ACont(:,1)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Ux [km/s^2]');
    subplot(3,1,2);
    plot(timeVec,ACont(:,2)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uy [km/s^2]');
    subplot(3,1,3);
    plot(timeVec,ACont(:,3)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('Uz [km/s^2]');
    xlabel('Time [days]');
    suptitle('Control acceleration in cart coordinates');
    set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');

    % Plot the modulus of the control acceleration
    figure;
    plot(timeVec,sqrt(ACont(:,1).^2 + ACont(:,2).^2 + ACont(:,3).^2)*AU2km/T_sid^2,'b','LineWidth',2.0);
    ylabel('U [km/s^2]');
    xlabel('Time [days]');
    suptitle('Modulus of control acceleration from Pseudo Equinoctial Shaping');
    
    %--------------------------------------------------------------------------
    % Pseudo elements plotting
    %--------------------------------------------------------------------------
    % First three
    figure;
    subplot(3,1,1);
    plot(LVec,Pseudo_Matrix(:,1),'b','linewidth',2);
    hold on;
    plot(LVec(1),Pseudo_Dep(1),'>',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    plot(LVec(end),Pseudo_Arr(1),'<',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    ylabel('p');
    subplot(3,1,2);
    plot(LVec,Pseudo_Matrix(:,2),'b','linewidth',2);
    hold on;
    plot(LVec(1),Pseudo_Dep(2),'>',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    plot(LVec(end),Pseudo_Arr(2),'<',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    ylabel('f');
    subplot(3,1,3);
    plot(LVec,Pseudo_Matrix(:,3),'b','linewidth',2);
    hold on;
    plot(LVec(1),Pseudo_Dep(3),'>',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    plot(LVec(end),Pseudo_Arr(3),'<',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    ylabel('g');
    suptitle('Pseudo elem from shaping');
    xlabel('LVec');
    
    % Last three
    figure;
    subplot(3,1,1);
    plot(LVec,Pseudo_Matrix(:,4),'b','linewidth',2);
    hold on;
    plot(LVec(1),Pseudo_Dep(4),'>',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    plot(LVec(end),Pseudo_Arr(4),'<',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    ylabel('h');
    subplot(3,1,2);
    plot(LVec,Pseudo_Matrix(:,5),'b','linewidth',2);
    hold on;
    plot(LVec(1),Pseudo_Dep(5),'>',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    plot(LVec(end),Pseudo_Arr(5),'<',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    ylabel('k');
    suptitle('Pseudo elem from shaping');
    xlabel('LVec');
    
end