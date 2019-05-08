function  plotTrajectory(DepBody,ArrBody,t_dep,TOF,thetaVec,R,Phi,Ucart)
% Function for plotting the low-thrust trajectory given the Departure and
% Arrival bodies parameters and the position vector in spherical
% coordinates of the spacecraft

% INPUT
% DepBody, ArrBody: departure and arrival bodies
% t_dep: departure time
% TOF: time of flight
% thetaVec: time vector of size (1 x timeStep)
% R: modulus of radius of the trajectory
% Phi: angle Phi of the trajectory
% Ucart: control acceleration in cartesian frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

%---------------------------------------------------------------------
% Plot thye orbits of Arrival and Departure bodies
%---------------------------------------------------------------------
% Angle for plot
theta_graph = [0:pi/100:2*pi];

% Departure body
p = DepBody.a*(1-DepBody.e^2);
r_graph = p ./ ( 1 + DepBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_geD] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

% Plot departure body
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

% Plot departure body
plot3(r_geA(1,:), r_geA(2,:), r_geA(3,:),'r','LineWidth',2.0);

%---------------------------------------------------------------------
% Find the low thrust trajectory and plot it
%---------------------------------------------------------------------
Xt = R.*cos(Phi).*cos(thetaVec);
Yt = R.*cos(Phi).*sin(thetaVec);
Zt = R.*sin(Phi);

% Plot the transfer
plot3(Xt,Yt,Zt,'k','LineWidth',2.0);
plot3(Xt(1),Yt(1),Zt(1),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','w',...
    'MarkerSize',8);
plot3(Xt(end),Yt(end),Zt(end),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',8);

%---------------------------------------------------------------------
% Plot position of Departure body at departure and Arrival body at arrival
%---------------------------------------------------------------------
% Find position of Departure Body at departure
M_dep = DepBody.M0 + DepBody.n*(t_dep - 0);
Theta_dep = M2theta(M_dep,DepBody.e);

% Departure body
p = DepBody.a*(1-DepBody.e^2);
r_dep = p ./ ( 1 + DepBody.e * cos(Theta_dep) );
X = cos(Theta_dep) .* r_dep;
Y = sin(Theta_dep) .* r_dep;
Z = 0 * r_dep;
r_pf = [X;Y;Z];
[r_ge] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

% Plot departure body
plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);

% Find position of Arrival Body at arrival
M_arr = ArrBody.M0 + ArrBody.n*(t_dep + TOF - 0);
Theta_arr = M2theta(M_arr,ArrBody.e);

% Arrival body
p = ArrBody.a*(1-ArrBody.e^2);
r_arr = p ./ ( 1 + ArrBody.e * cos(Theta_arr) );
X = cos(Theta_arr) .* r_arr;
Y = sin(Theta_arr) .* r_arr;
Z = 0 * r_arr;
r_pf = [X;Y;Z];
[r_ge] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);

% Plot arrival body
plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a cube
% side = 0.05;
% vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1]*side;
% % Add the x,y,z coordinates of the point in the graph where to plot the
% % cube
% vert(:,1) = vert(:,1) + Xt(35) - side/2;
% vert(:,2) = vert(:,2) + Yt(35) - side/2;
% vert(:,3) = vert(:,3) + Zt(35) - side/2;
% fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
% patch('Vertices',vert,'Faces',fac,...
% 'FaceVertexCData',hsv(8),'FaceColor','interp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
% Legend,axis,title
%---------------------------------------------------------------------
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

% Plot acceleration vector
molt = 1e5;
Ucart = molt*Ucart;
quiver3(Xt,Yt,Zt,Ucart(1,:),Ucart(2,:),Ucart(3,:),'b','LineWidth',2.0);

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
 
