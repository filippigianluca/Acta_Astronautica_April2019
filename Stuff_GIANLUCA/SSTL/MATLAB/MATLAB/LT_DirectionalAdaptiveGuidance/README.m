% Directional Adaptive Guidance
% Reference: Ruggiero, Pergola, Marcuccio, Andrenucci, "Low-Thrust
% Manuevers for the Efficient Correction of Orbital Elements", 32nd
% International Electric Propulsion Conference, 2011

% Input: initial orbital elements, final orbital elements, tolerance for
% the efficiency of the variation of the orbital elements
% Output: control history, DeltaV, ToF

% 
% - This is a targeting algorithm, it does not guarantees optimality
% - The thrust vector is a linear combination of LT accelerations that give
% the maximum instantanoues variation of each orbital element 
% - When the manuever is not efficient, the engine is switched off