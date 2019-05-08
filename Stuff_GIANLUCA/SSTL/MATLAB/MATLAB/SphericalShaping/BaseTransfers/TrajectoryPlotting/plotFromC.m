%% File to print results obtained from C++ simulations
%% File for single trajectory and for pork chop plots for the 4 test cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data and parameters
% clc;
% clear all;
% close all;

% Add folders containing subfunctions
%addpath(genpath('SphericalShapingFunctions'));
addpath('SphericalShapingFunctions');
%addpath(genpath('TrajectoryPlotting'));
addpath('TrajectoryPlotting');

% Load the data
LoadData;

%% Select bodies

% Define Earth as departure body
DepBody = body(1);

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

%% Load results
R = csvread('SphericalShapingR.csv');
Phi = csvread('SphericalShapingPhi.csv');
thetaVec =  csvread('SphericalShapingthetaVec.csv');
u_mat =  csvread('SphericalShapingu_mat.csv');

%% Function for plotting the low-thrust trajectory given the Departure and
% Arrival bodies parameters and the position vector in spherical
% coordinates of the spacecraft
%------------------------------------------------------------------------%

% Angle for plot
theta_graph = [0:pi/200:2*pi];

% Low-thrust trajectory
Xt = R.*cos(Phi).*cos(thetaVec);
Yt = R.*cos(Phi).*sin(thetaVec);
Zt = R.*sin(Phi);

% Departure body
p = DepBody.a*(1-DepBody.e^2);
r_graph = p ./ ( 1 + DepBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geD] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

figure;
plot3(r_geD(1,:), r_geD(2,:), r_geD(3,:),'g','LineWidth',2.0);
hold on;
grid on;


% Arrival body
p = ArrBody.a*(1-ArrBody.e^2);
r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geA] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);

plot3(r_geA(1,:), r_geA(2,:), r_geA(3,:),'r','LineWidth',2.0);

% Plot the transfer
plot3(Xt,Yt,Zt,'k','LineWidth',2.0);
plot3(Xt(1),Yt(1),Zt(1),'gx','LineWidth',3.0);
plot3(Xt(end),Yt(end),Zt(end),'rx','LineWidth',3.0);

% Legend
leg1 = strcat('Departure orbit: ',DepBody.name);
leg2 = strcat('Arrival orbit: ',ArrBody.name);
legend(leg1,leg2,'Transfer trajectory','location','southoutside');

%Plot axes
if ArrBody.a > 5
    lim = 35;
else
    lim = 4;
end
line([0 lim], [0 0], [0 0]); text(lim, 0, 0, 'X');
line( [0 0], [0 lim], [0 0]); text( 0, lim, 0, 'Y');
line( [0 0], [0 0], [0 lim]); text( 0, 0, lim, 'Z');

%Plot axis
axis equal;
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
tit = strcat('Low thrust transfer trajectory: ',DepBody.name, '-',ArrBody.name);
title(tit);
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
view(3);


%% Plot the control as a function of angle theta (not time!!)
figure;
subplot(3,1,1);
plot(thetaVec*r2d,u_mat(1,:),'b','LineWidth',2.0);
ylabel('Ut [km/s^2]');
subplot(3,1,2);
plot(thetaVec*r2d,u_mat(2,:),'b','LineWidth',2.0);
ylabel('Un [km/s^2]');
subplot(3,1,3);
plot(thetaVec*r2d,u_mat(3,:),'b','LineWidth',2.0);
ylabel('Uh [km/s^2]');
xlabel('Theta [deg]');
suptitle('Control acceleration');
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');


%% Pork chop plot for the Grid Search
DepartureDate = csvread('SphericalShapingDepartureDates.csv');
TransferTime = csvread('SphericalShapingTransferTimes.csv');
DV_matrix =  csvread('SphericalShapingDeltaV.csv');
X = [];
Y = [];

%Plot of the data
for i=1:length(TransferTime)
   X(i,:) = DepartureDate'; 
end
for i=1:length(DepartureDate)
   Y(:,i) = TransferTime; 
end
figure;
contourf(X,Y,DV_matrix*AU2km/T_sid);
c = colorbar;
title('DeltaV obtained [km/s]');
xlabel('Departure date [MJD2000]');
ylabel('Transfer time [days]');
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');

% Filter out DV higher than 50 km/s
for i=1:length(TransferTime)
    for j=1:length(DepartureDate)
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


%% Plot movie of trajectory

fH = figure;
figure(fH);
hold on;
grid on;

% Set axes limit
if ArrBody.a > 5
    lim = 35;
else
    lim = 4;
end

%Mov = getframe;
for i=1:length(thetaVec)
    figure(fH)
    clf, hold on, view(3), grid, camproj('perspective')
    line([0 lim], [0 0], [0 0]); text(lim, 0, 0, 'X');
    line( [0 0], [0 lim], [0 0]); text( 0, lim, 0, 'Y');
    line( [0 0], [0 0], [0 lim]); text( 0, 0, lim, 'Z');
    plot3(r_geD(1,:), r_geD(2,:), r_geD(3,:),'g','LineWidth',2.0);
    plot3(r_geA(1,:), r_geA(2,:), r_geA(3,:),'r','LineWidth',2.0);
    plot3(Xt,Yt,Zt,'k--','LineWidth',1.0); 
    plot3(Xt(1),Yt(1),Zt(1),'gx','LineWidth',3.0);
    plot3(Xt(end),Yt(end),Zt(end),'rx','LineWidth',3.0);
    plot3(Xt(1:i),Yt(1:i),Zt(1:i),'k','LineWidth',2.0); 
    plot3(Xt(i),Yt(i),Zt(i),'kx','LineWidth',3.0); 
    drawnow;
    Mov(i) = getframe(fH);
end

movie(Mov,1,60)

%Plot axis
axis equal;
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
tit = strcat('Low thrust transfer trajectory: ',DepBody.name, '-',ArrBody.name);
title(tit);
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
view(3);