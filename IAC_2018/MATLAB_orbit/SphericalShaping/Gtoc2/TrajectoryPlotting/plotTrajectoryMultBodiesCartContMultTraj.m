function  plotTrajectoryMultBodiesCartContMultTraj(DepBody,ArrivalBody,depDate,arrDate,R,Ucart)
% Function for plotting the low-thrust trajectory given the Departure and
% Arrival bodies parameters and the position vector in spherical
% coordinates of the spacecraft, for the GTOC2 asteroids.
% The legs that have been computed successfully by the spherical shaping
% are plotted (even if not all the legs have been computed)
% Modified for plotting also the control acceleration
% 
% Version for tha different segments contained in different layers of the
% matrix (multidimensional matrix)
%
% INPUT: 
% DepBody,ArrivalBody: structures containing the parameters on the
% bodies
% depDate, arrDate: vector with departure and arrival
% R,Ucart: matrices containing the values of R,Ucart in cartesian coord
% for each shaped leg (if in the middle of two shaped legs there is a leg
% which has not been shaped, then the value is 0 and this cause no problem)
%
% OUTPUT:
% no output, just graph plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the orbits

%---------------------------------------------------------------------
% Plot all the orbits for the bodies as dashed lines
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
[r_ge] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

% Plot departure body
figure;
plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'g--','LineWidth',2.0);
hold on;
grid on;

% Arrival bodies
for i=1:4
    ArrBody = ArrivalBody(i);
    p = ArrBody.a*(1-ArrBody.e^2);
    r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
    X = cos(theta_graph) .* r_graph;
    Y = sin(theta_graph) .* r_graph;
    Z = 0 * r_graph;
    r_pf = [X;Y;Z];
    [r_ge] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);
    
    % Plot Arrival bodies
    plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'r--','LineWidth',1.0);
end

%---------------------------------------------------------------------
% Plot the low-thrust trajectory
%---------------------------------------------------------------------
for i=1:4
    
    % Extract trajectory in cartesian coordinates
    Xt(:) = R(i,1,:);
    Yt(:) = R(i,2,:);
    Zt(:) = R(i,3,:);
    
    % Plot the transfer
    plot3(Xt,Yt,Zt,'k','LineWidth',2.0);
    plot3(Xt(1),Yt(1),Zt(1),'o',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','w',...
        'MarkerSize',10);
    plot3(Xt(end),Yt(end),Zt(end),'o',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10);
    
    % Plot acceleration vector
    molt = 1e5;
    Ucart = molt*Ucart;
    Ucartx(:) = Ucart(i,1,:);
    Ucarty(:) = Ucart(i,2,:);
    Ucartz(:) = Ucart(i,3,:);
    quiver3(Xt,Yt,Zt,Ucartx,Ucarty,Ucartz,'b','LineWidth',2.0);

end

%---------------------------------------------------------------------
% Plot the orbits of the body with departure and arrival times
%---------------------------------------------------------------------
% Find position and velocity of Departure Body at initial departure
M_dep = DepBody.M0 + DepBody.n*(depDate(1)-DepBody.t0);
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

% Arrival bodies
for i=1:3
    
    % Find arrival body of the leg
    ArrBody = ArrivalBody(i);
    
    % Arrival
    M_arr = ArrBody.M0 + ArrBody.n*(arrDate(i)-ArrBody.t0);
    Theta_arr = M2theta(M_arr,ArrBody.e);
    
    % Departure
    M_dep = ArrBody.M0 + ArrBody.n*(depDate(i+1)-ArrBody.t0);
    Theta_dep = M2theta(M_dep,ArrBody.e);
    
    while (Theta_dep < Theta_arr)
        Theta_dep = Theta_dep + 2*pi;
    end
    
    % Angle for plotting
    theta_graph = [Theta_arr:pi/100:Theta_dep];
    
    p = ArrBody.a*(1-ArrBody.e^2);
    r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
    X = cos(theta_graph) .* r_graph;
    Y = sin(theta_graph) .* r_graph;
    Z = 0 * r_graph;
    r_pf = [X;Y;Z];
    [r_ge] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);
    
    % Plot Arrival bodies
    plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'r','LineWidth',2.0);
    plot3(r_ge(1,1), r_ge(2,1), r_ge(3,1),'o',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',8);
    plot3(r_ge(1,end), r_ge(2,end), r_ge(3,end),'o',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',8);
    
end

% For last body
% Find arrival body of the leg
ArrBody = ArrivalBody(4);

% Arrival
M_arr = ArrBody.M0 + ArrBody.n*(arrDate(4)-ArrBody.t0);
Theta_arr = M2theta(M_arr,ArrBody.e);

% Parameters of last body
p = ArrBody.a*(1-ArrBody.e^2);
r_arr = p ./ ( 1 + ArrBody.e * cos(Theta_arr) );
X = cos(Theta_arr) .* r_arr;
Y = sin(Theta_arr) .* r_arr;
Z = 0 * r_arr;
r_pf = [X;Y;Z];
[r_ge] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);

% Plot last arrival body
plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);

%---------------------------------------------------------------------
% Legend,axis,title
%---------------------------------------------------------------------
% Legend
leg1 = strcat('Departure orbit: ',DepBody.name);
legend(leg1,'location','southoutside');

%Plot axes
lim = 7;
line([0 lim], [0 0], [0 0]); text(lim, 0, 0, 'X');
line( [0 0], [0 lim], [0 0]); text( 0, lim, 0, 'Y');
line( [0 0], [0 0], [0 lim]); text( 0, 0, lim, 'Z');

%Plot axis
axis equal;
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
tit = strcat('Trajectory ',DepBody.name, '-',ArrivalBody(1).name, '-',ArrivalBody(2).name, '-',ArrivalBody(3).name, '-',ArrivalBody(4).name);
title(tit);
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
view(3);

end
 
