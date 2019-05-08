function plot_planet_orbits(t0,TOF,s,varargin)
% Plots the orbits of the planets (or asteroids) in the solar system in
% the given interval of time [t0, t0+TOF], and optionally plots also one or
% more trajectories, given point by point.
%
%   plot_planet_orbits(t0, TOF, s [, r1 [, r2 [, ...]]])
%
% Plots in a new figure, in AU.
%
% INPUT
%   t0      Initial time of integration [d, MJD2000].
%   TOF     Time of integration [d].
%   s       Vector containing the id of all the bodies (planets and/or
%           asteroids) to be plotted.
%   r1, ... Matrix containing other trajectories to be plotted [km]. Each
%           trajectory is given in ecliptic cartesian coordinates, one
%           point of the trajectory for each row.
%
% OUTPUT
%   (none)
%
% FUNCTIONS CALLED
%   two_body_dynamics, EphSS
%
% Matteo Ceriotti, 06-02-2007

mu=astro_constants(4);
% Initialises figure
figure;
hold on;
line(0,0,0,'Marker','*','Color','r');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
AU=astro_constants(2);
% Initialise integrator used to plot trajectories
ode45options = odeset('abstol',1e-11,'reltol',1e-9);

% Plots all the planets from time to time+tof
for i=1:length(s)
    [xpi,vpi]=EphSS(s(i),t0); % Position and velocity of the arrival planet
    [ti,yi]=ode45(@two_body_dynamics,[0, TOF*86400],[xpi,vpi],ode45options,mu);
    line(yi(:,1)/AU,yi(:,2)/AU,yi(:,3)/AU,'Color',[s(i)/max(s) 0.5 0]);
    line(yi(1,1)/AU,yi(1,2)/AU,yi(1,3)/AU,'Color',[s(i)/max(s) 0.5 0],'Marker','.');
    line(yi(end,1)/AU,yi(end,2)/AU,yi(end,3)/AU,'Color',[s(i)/max(s) 0.5 0],'Marker','.');
end

for i=1:length(varargin)
    line(varargin{i}(:,1)/AU,varargin{i}(:,2)/AU,varargin{i}(:,3)/AU,'Color','b');
    line(varargin{i}(1,1)/AU,varargin{i}(1,2)/AU,varargin{i}(1,3)/AU,'Marker','*') % Marker at departure point
    line(varargin{i}(end,1)/AU,varargin{i}(end,2)/AU,varargin{i}(end,3)/AU,'Marker','o') % Marker at arrival point
end

axis equal;

return