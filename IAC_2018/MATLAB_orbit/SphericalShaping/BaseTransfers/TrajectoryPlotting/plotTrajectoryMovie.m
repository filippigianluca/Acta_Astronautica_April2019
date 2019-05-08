function  plotTrajectoryMovie(DepBody,ArrBody,thetaVec,R,Phi,mu,t_dep,TOF,Ucart)
% Function for plotting the low-thrust trajectory given the Departure and
% Arrival bodies parameters and the position vector in spherical
% coordinates of the spacecraft
% This function plots the movie of the trajectory

% INPUT
% DepBody, ArrBody: departure and arrival bodies
% t_dep: departure time
% TOF: time of flight
% thetaVec: time vector of size (1 x timeStep)
% R: modulus of radius of the trajectory
% Phi: angle Phi of the trajectory
% Ucart: control acceleration in cartesian frame
% mu: gravitational parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find angles theta in perifocal frame for initial and final state of Departure and Arrival bodies

M_dep = DepBody.M0 + DepBody.n*(t_dep-0);
Theta_dep0 = M2theta(M_dep,DepBody.e);

M_arr = DepBody.M0 + DepBody.n*(t_dep+TOF-0);
Theta_depf = M2theta(M_arr,DepBody.e);

if (Theta_depf < Theta_dep0)
    Theta_depf = Theta_depf + 2*pi;
end
DepP = 2*pi*sqrt(DepBody.a^3/mu);
nr = floor(TOF/DepP);
Theta_depf = Theta_depf + 2*pi*nr;

M_dep = ArrBody.M0 + ArrBody.n*(t_dep-0);
Theta_arr0 = M2theta(M_dep,ArrBody.e);

M_arr = ArrBody.M0 + ArrBody.n*(t_dep+TOF-0);
Theta_arrf = M2theta(M_arr,ArrBody.e);

if (Theta_arrf < Theta_arr0)
    Theta_arrf = Theta_arrf + 2*pi;
end
ArrP = 2*pi*sqrt(ArrBody.a^3/mu);
nr = floor(TOF/ArrP);
Theta_arrf = Theta_arrf + 2*pi*nr;

%% Set vectors with Departure and Arrival bodies states from [0-2 pi]

% Angle for plot
theta_graph = [0:pi/200:2*pi];

% Low-thrust trajectory
Xt = R.*cos(Phi).*cos(thetaVec);
Yt = R.*cos(Phi).*sin(thetaVec);
Zt = R.*sin(Phi);
l = length(thetaVec);

% Departure body
p = DepBody.a*(1-DepBody.e^2);
r_graph = p ./ ( 1 + DepBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geD] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

% Arrival body
p = ArrBody.a*(1-ArrBody.e^2);
r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geA] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);

%% Set vectors with Departure and Arrival bodies states from [theta_arr - theta_dep]

% Departure body
theta_graph = linspace(Theta_dep0,Theta_depf,l);

p = DepBody.a*(1-DepBody.e^2);
r_graph = p ./ ( 1 + DepBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geDt] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

% Arrival body
theta_graph = linspace(Theta_arr0,Theta_arrf,l);

p = ArrBody.a*(1-ArrBody.e^2);
r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geAt] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);

%% Plot movie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For cube plotting
% side = 0.1;
% vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1]*side;
% fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set handle
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

% Multiplicator for represenation of acceleration
molt = 1e5;

% Plot at each theta step and save for movie
for i=1:length(thetaVec)
    figure(fH)
    clf, hold on, view([1 1 1]), grid, camproj('perspective')
    line([0 lim], [0 0], [0 0]); text(lim, 0, 0, 'X');
    line( [0 0], [0 lim], [0 0]); text( 0, lim, 0, 'Y');
    line( [0 0], [0 0], [0 lim]); text( 0, 0, lim, 'Z');
    plot3(r_geD(1,:), r_geD(2,:), r_geD(3,:),'g--','LineWidth',1.0);
    plot3(r_geA(1,:), r_geA(2,:), r_geA(3,:),'r--','LineWidth',1.0);
    plot3(Xt,Yt,Zt,'k--','LineWidth',1.0); 
%     plot3(Xt(1),Yt(1),Zt(1),'gx','LineWidth',3.0);
%     plot3(Xt(end),Yt(end),Zt(end),'rx','LineWidth',3.0);
    plot3(r_geDt(1,1:i), r_geDt(2,1:i), r_geDt(3,1:i),'g','LineWidth',2.0);
    plot3(r_geAt(1,1:i), r_geAt(2,1:i), r_geAt(3,1:i),'r','LineWidth',2.0);
    plot3(r_geDt(1,i), r_geDt(2,i), r_geDt(3,i),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);
    plot3(r_geAt(1,i), r_geAt(2,i), r_geAt(3,i),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);
    plot3(Xt(1:i),Yt(1:i),Zt(1:i),'k','LineWidth',2.0); 
    plot3(Xt(i),Yt(i),Zt(i),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',8);
%   quiver3(Xt(1:i),Yt(1:i),Zt(1:i),Ucart(1,(1:i)),Ucart(2,(1:i)),Ucart(3,(1:i)),'b','LineWidth',1.0);
    quiver3(Xt(i),Yt(i),Zt(i),molt*Ucart(1,(i)),molt*Ucart(2,(i)),molt*Ucart(3,(i)),'b','LineWidth',2.0);
    

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a cube
% Add the x,y,z coordinates of the point in the graph where to plot the
% cube
% vert(:,1) = vert(:,1) + Xt(i) - side/2;
% vert(:,2) = vert(:,2) + Yt(i) - side/2;
% vert(:,3) = vert(:,3) + Zt(i) - side/2;
% patch('Vertices',vert,'Faces',fac,...
% 'FaceVertexCData',hsv(6),'FaceColor','flat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    drawnow;
    Mov(i) = getframe(fH);
end

movie(Mov,1,4)


%Plot axis
axis equal;
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
tit = strcat('Low thrust transfer trajectory: ',DepBody.name, '-',ArrBody.name);
title(tit);
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
view([1 1 1]);
end

