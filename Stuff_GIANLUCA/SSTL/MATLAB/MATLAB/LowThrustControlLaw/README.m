% Low-thrust laws for the variation of specific orbital
% elements using low-thrust propulsion.


%% Compute final orbital elements and DeltaV from initial orbital elements and ToF

% Input: initial orbital elements, ToF
% Output: variation of each orbital element in that ToF and DeltaV

% Possible application: pruning of the search space based on variation of orbital
% elements that is possible to obtain in a given ToF

% File "OE0ToF_2_OEfDeltaV.m"


%% Compute DeltaV and ToF from initial and final orbital elements

% Input: initial orbital elements, final orbital elements
% Output: DeltaV, ToF

% Possible application: low-fidelity model for the evaluation of the DeltaV
% of orbital transfer

% File "OE0OEf_2_DeltaVToF.m"


%% Files with plot and numerical comparison for each law

% -------------------------------------------------------------------------
% SEMIMAJOR AXIS
% Folder "SemimajorAxis"
% -------------------------------------------------------------------------
% 1. Maximum rate of change of semimajor axis: thrust in the tangential
% direction (flight path angle):
% a_Optimal.m
%
% 2. Change of semiamjor axis and inclination according to Edelbaum:
% Edelbaum_a_i_main.m
%
% 3. Change of semimajor axis and inclination with Edelbaum's law and with
% constraint on the maximum value of the semimajor axis during the transfer
% (Kechichian):
% Kechichian_a_i_constrained.m
%
% 4. Minimum time variation of semimajor axis and right ascension:
% Kechichian_a_RAAN_main.m
%
% 5. Variation of semimajor axis, eccentricity and argument of the perigee:
% Fernandes_main.m
%
% 6. Variation of semimajor axis with no variation of eccentricity:
% Burt_Delta_a_NoDelta_e.m
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ECCENTRICITY
% Folder "Eccentricity"
% -------------------------------------------------------------------------
% 1. Maximum rate of change of the eccentricity:
% e_Optimal.m
%
% 2. Variation of eccentricity with no variation of semimajor axis:
% Burt_Delta_e_NoDelta_a.m
% 
% 3. Simultanous variation of eccentricity and inclination
% Pollard_main_e_i.m
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% INCLINATION
% Folder "Inclination"
% -------------------------------------------------------------------------
% 1. Maximum variation of inclination:
% i_Optimal.m
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ARGUMENT OF THE PERIGEE
% Folder "ArgumentOfPerigee"
% ------------------------------------------------------------------------- 
% 1. Variation of the argument of perigee without variation of semimarjor
% axis and eccentricity:
% Pollard_main_omega.m
%
% 2. Maximum variation of the argument of perigee according to Petropolous:
% omega_Optimal_Petropolous.m
%
% 3. Maximum variation of the argument of perigee according to Ruggiero:
% omega_Optimal_Ruggiero.m
%
% 4. Variation of the argument of the perigee without variation of
% semimajor axis and eccentricity, according to Burt (1):
% Burt_Delta_omega_NoDelta_a_e.m
%
% 5. Variation of the argument of the perigee without variation of
% semimajor axis and eccentricity, Burt (2):
% Burt_Delta_omega_NoDelta_a_e_2.m
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% RIGHT ASCENSION OF THE ASCENDING NODE
% Folder "RightAscensionAscendingNode"
% -------------------------------------------------------------------------
% 1. Maximum variation of right ascension:
% RAAN_Optimal.m
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
