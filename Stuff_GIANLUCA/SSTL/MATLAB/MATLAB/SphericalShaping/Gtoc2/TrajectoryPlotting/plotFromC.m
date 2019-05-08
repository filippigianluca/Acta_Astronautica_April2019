%% File to print results obtained from C++ simulations
%% File for the GTOC2 asteroids trajectory
%% Contrary to the correspondent plot for the Matlab file, here the trajectory is plotted only if all the 4 legs have been shaped
%% (if file from C++ are provided, which means that a trajectory has been found)
%%-------------------------------------------------------------------%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and parameters
% clc;
% clear all;
% close all;

% Define conversion parameters
d2r = pi/180;
r2d = 180/pi;
AU2km = 149597870.66;                         % [km]
km2AU = 1/AU2km;

% Set constants
T_sid = 86164.1004;                           % [s]
mu_S = 1.327178*1e11*km2AU^3*T_sid^2;         % [AU^3/day^2]

% Load the asteroids data
AstComb;

% Load file with asteroids data (ATTENTION: ANGLE IN DEGREES)
[data,f] = readFileImportdata('ast_ephem.txt');

%% Select bodies

% Define Earth as departure body
Earth.name = 'Earth';
Earth.a = 0.999988049532578;
Earth.e = 1.671681163160e-02;
Earth.i = 0.8854353079654e-03*d2r;
Earth.Omega = 175.40647696473*d2r;
Earth.w = 287.61577546182*d2r;
Earth.M0 = 257.60683707535*d2r;
Earth.n = sqrt(mu_S/Earth.a^3);     
Earth.t0 = 54000;

DepBody = Earth;

choice = menu('Select the asteroids combination','1','2','3','4','5','6','7','8','9','10','11','12');
switch choice
    case 1
       comb = 1; 
    case 2
       comb = 2; 
    case 3
       comb = 3; 
    case 4
       comb = 4; 
    case 5
       comb = 5; 
    case 6
       comb = 6; 
    case 7
       comb = 7; 
    case 8
       comb = 8; 
    case 9
       comb = 9; 
    case 10
       comb = 10; 
    case 11
       comb = 11; 
    case 12
       comb = 12; 
end

AST = ASTEROIDS(comb,:);
depDate_vec = depDate_vector(comb,:);
arrDate_vec = arrDate_vector(comb,:);

% Search for the desired asteroids in the database and insert in the
% structure ArrBody
for i=1:4
    for j=1:max(size(data))-1
        if data(j,1)==AST(i)
            ArrivalBody(i).name = num2str(AST(i));
            ArrivalBody(i).a = data(j,2);
            ArrivalBody(i).e = data(j,3);
            ArrivalBody(i).i = data(j,4)*d2r;
            ArrivalBody(i).Omega = data(j,5)*d2r;
            ArrivalBody(i).w = data(j,6)*d2r;
            ArrivalBody(i).M0 = data(j,7)*d2r;
            ArrivalBody(i).n = sqrt(mu_S/ArrivalBody(i).a^3);
            ArrivalBody(i).t0 = data(j,8);
        end
    end
end

%% Load files

str = '';
uT = textread('uT.dat', str);
uN = textread('uN.dat', str);
uH = textread('uH.dat', str);
thetalegs = textread('thetalegs.dat', str);
Rlegs = textread('Rlegs.dat', str);
Philegs = textread('Philegs.dat', str);

% Add 0 at the end of each row so that following algorithm doesnt crash
% in correspondence of the line with all nonzero elements
uT(:,end+1) = 0;
uN(:,end+1) = 0;
uH(:,end+1) = 0;
thetalegs(:,end+1) = 0;
Rlegs(:,end+1) = 0;
Philegs(:,end+1) = 0;

%break
%% Function for plotting the low-thrust trajectory given the Departure and
% Arrival bodies parameters and the position vector in spherical
% coordinates of the spacecraft
%------------------------------------------------------------------------%

% Angle for plot
theta_graph = [0:pi/200:2*pi];

% Departure body
p = DepBody.a*(1-DepBody.e^2);
r_graph = p ./ ( 1 + DepBody.e * cos(theta_graph) );
X = cos(theta_graph) .* r_graph;
Y = sin(theta_graph) .* r_graph;
Z = 0 * r_graph;
r_pf = [X;Y;Z];
[r_ge] = KeplElem2rvMatrix(DepBody.Omega, DepBody.i, DepBody.w, r_pf);

figure;
plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'g','LineWidth',1.0);
hold on;
grid on;


% Arrival body
for i=1:4
    ArrBody = ArrivalBody(i);
    p = ArrBody.a*(1-ArrBody.e^2);
    r_graph = p ./ ( 1 + ArrBody.e * cos(theta_graph) );
    X = cos(theta_graph) .* r_graph;
    Y = sin(theta_graph) .* r_graph;
    Z = 0 * r_graph;
    r_pf = [X;Y;Z];
    [r_ge] = KeplElem2rvMatrix(ArrBody.Omega, ArrBody.i, ArrBody.w, r_pf);
    
    plot3(r_ge(1,:), r_ge(2,:), r_ge(3,:),'r','LineWidth',1.0);
end

for i=1:4
    % Select leg (excluding first elements which is only leg indicator)
    thetaVec = thetalegs(i,2:end);
    R = Rlegs(i,2:end);
    Phi = Philegs(i,2:end);
    
    % Find index correspondent to first zero element
    a = find(thetaVec==0);
    a = a(1) - 1;
    
    % Exclude all the zero elements from the vectors
    thetaVec = thetaVec(1:a);
    R = R(1:a);
    Phi = Phi(1:a);
    
    Xt = R.*cos(Phi).*cos(thetaVec);
    Yt = R.*cos(Phi).*sin(thetaVec);
    Zt = R.*sin(Phi);
    
    % Plot the transfer
    plot3(Xt,Yt,Zt,'k','LineWidth',2.0);
    plot3(Xt(1),Yt(1),Zt(1),'kx','LineWidth',3.0);
    plot3(Xt(end),Yt(end),Zt(end),'kx','LineWidth',3.0);
end

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

%% Plot the control of the 4 legs


for i=1:4
    % Select leg (excluding first elements which is only leg indicator)
    uTplot = uT(i,2:end);
    uNplot = uN(i,2:end);
    uHplot = uH(i,2:end);
    thetaVec = thetalegs(i,2:end);
    
    % Find index correspondent to first zero element
    a = find(uTplot==0);
    a = a(1) - 1;
    
    % Exclude all the zero elements from the vectors
    uTplot = uTplot(1:a);
    uNplot = uNplot(1:a);
    uHplot = uHplot(1:a);
    thetaVec = thetaVec(1:a);
    
    figure;
    subplot(3,1,1);
    plot(thetaVec*r2d,uTplot,'b','LineWidth',2.0);
    ylabel('Ut [km/s^2]');
    subplot(3,1,2);
    plot(thetaVec*r2d,uNplot,'b','LineWidth',2.0);
    ylabel('Un [km/s^2]');
    subplot(3,1,3);
    plot(thetaVec*r2d,uHplot,'b','LineWidth',2.0);
    ylabel('Uh [km/s^2]');
    xlabel('Theta [deg]');
    tit = strcat('Control acceleration for leg N', num2str(i));
    suptitle(tit);
    set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
end


