% Control laws for the variation of semimajor axis, inclination and right
% ascension of the ascending node of circular orbit in a given time of
% flight under the effect of J2.
% Developed for the deployment of satellites from the injection orbit of
% the launcher to the operational orbit of a MEO constellation.
%
% Reference: M. Di Carlo, L. Ricciardi, M. Vasile, "Multi-objective
% optimisation of constellation deployment using low-thrust propulsion",
% SPACE 2016, Long Beach, CA
%
% Folders:
%
% -------------------------------------------------------------------------
% FirstStrategy
% -------------------------------------------------------------------------
% Variation of RAAN at low altitude to exploit J2, then variation of
% semimajor axis and inclination using propulsion on two thrust arc per
% revolution

% -------------------------------------------------------------------------
% SecondStrategy
% -------------------------------------------------------------------------
% Variation semimajor axis and inclination using propulsion on two thrust
% arc per revolution (no change of RAAN during this phase due to
% propulsion), then variation of RAAN using 90 deg elevation thrust on two
% thrust arc per revolution (thrust arcs are placed in appropriate
% positions on orbit so that no variation of inclination takes place during
% this phase)

% -------------------------------------------------------------------------
% ThirdStrategy
% -------------------------------------------------------------------------
% Variation of semimajor axis using thrust with elevation angle equal to
% zero on two thrust arcs per revolution (RAAN changes due to J2 and
% because of the variation of semimajor axis), then
% change of thrust profile. In the second part the engine is on with
% elevation angle equal to 90 deg during two thrust arcs per revolution to
% change the inclination (RAAN changes due to J2 and because of the
% variation of inclination). Final values of RAAN matches final desired
% RAAN.

