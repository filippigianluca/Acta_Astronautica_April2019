%% ground_station_list: this function generates structures relative to the 
%                       elliptical model used to represent the Earth and to various ground stations' parameters.
%
%% Inputs:
% * () : no input
%
%
%% Output:
% * GS : structure containing ground stations' data
% * Ell_model : structure containing Earth ellipsoid's parameters
%
% Author: Francesco Torre
% email: francesco.torre@strath.ac.uk

function [GS, ell_model] = ground_station_list()

%% Ellipsoid model

ell_model(1).name   = 'WSG84';
ell_model(1).R      = 6378.137; % Km
ell_model(1).f      = 1/298.257223563;

ell_model(2).name   = 'GEM-10B';
ell_model(2).R      = 6378.138; % Km
ell_model(2).f      = 1/298.257;



%% Ground stations

% Ground station 1 - DSA-1, New Norcia, ESA
GS(1).name      = 'DSA-1';
GS(1).phi       = -031  +02/60  +53.61/3600;    % deg - latitude
GS(1).lambda    = +116  +11/60  +29.40/3600;    % deg - longitude
GS(1).h         = 0.25226;   % Km - altitude/WSG84
GS(1).E0        = 3;  %check!     % deg - minimum elevation
GS(1).marker    = '+';
GS(1).colour    = 'r';
GS(1).Rn_err    = 10; %check!     % Range error [Km]
GS(1).Az_err    = 6e-3; %check!   % Azimuth error [degrees]
GS(1).El_err    = 6e-3; %check!   % Elevation error [degrees]
GS(1).Rr_err    = 1e-2; %check!   % Range rate error [Km/s]

% Ground station 2 - DSA-2, Cebreros, ESA
GS(2).name      = 'DSA-2';
GS(2).phi       = +040  +27/60  +09.68/3600;    % deg - latitude
GS(2).lambda    = +004  +22/60  +03.18/3600;    % deg - longitude
GS(2).h         = 0.79410;   % Km - altitude/WSG84
GS(2).E0        = 3;  %check!    % deg - minimum elevation
GS(2).marker    = 'o';
GS(2).colour    = 'g';
GS(2).Rn_err    = 10; %check!     % Range error [Km]
GS(2).Az_err    = 6e-3; %check!   % Azimuth error [degrees]
GS(2).El_err    = 6e-3; %check!   % Elevation error [degrees]
GS(2).Rr_err    = 1e-3; %check!   % Range rate error [Km/s]

% Ground station 3 - DSA-3, Malargue, ESA
GS(3).name      = 'DSA-3';
GS(3).phi       = +035  +46/60  +33.63/3600;    % deg - latitude
GS(3).lambda    = +069  +23/60  +53.51/3600;    % deg - longitude
GS(3).h         = 1.550;   % Km - altitude/WSG84
GS(3).E0        = 3;  %check!    % deg - minimum elevation
GS(3).marker    = 'x';
GS(3).colour    = 'b';
GS(3).Rn_err    = 10; %check!     % Range error [Km]
GS(3).Az_err    = 6e-3; %check!   % Azimuth error [degrees]
GS(3).El_err    = 6e-3; %check!   % Elevation error [degrees]
GS(3).Rr_err    = 1e-3; %check!   % Range rate error [Km/s]



