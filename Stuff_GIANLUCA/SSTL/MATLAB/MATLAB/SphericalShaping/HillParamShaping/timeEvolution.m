function [timeVec] = timeEvolution(thetaVec,TPrime,t_dep)

% timeEvolution
% Function that calculates the temporal story of the trajectory given in
% input the vectors of theta and time derivative respect to angle.
% Use the fact that the theta step is constant
%
% To decrease the error the derivative is approximated as deriv = (TPrime(i-1) + TPrime(i))/2
% This has shown to yield an higher accuracy than using just TPrime(i-1) or TPrime(i)
%
% INPUT
% thetaVec: vector with angle theta [rad]
% TPrime: vector with dT/dtheta     [days/rad]
% t_dep: departure day              [JD,MJD,MJD2000]
%
% OUTPUT
% timeVec: vector with time elements of the trajectory (the time
% correspondent at each theta point) [JD,MJD,MJD2000]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute temporal evolution

% Compute thetaStep and first time element
thetaStep = thetaVec(2)- thetaVec(1);

timeVec = [];
timeVec(1) = t_dep;

% Compute all the next time element with a finite derivative approximation
for i=2:length(thetaVec)
 deriv = (TPrime(i-1) + TPrime(i))/2;
 timeVec(i) = timeVec(i-1) + thetaStep*deriv;
end


end

