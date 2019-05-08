function [FunVal] = DirTranscFunctionFMINCONMatObj(y,N,tVec)

%DirTranscFunctionFMINCONMatObj
% Function that calculates discretized performance index for the direct
% transcription problem of the optimal control problem applied to low
% thrust trajectory optimization
%
% INPUT
% y: array comprising the state and the control at each moment of time
% y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
% a 3x1 vector (velocity) and u is a 3x1 vector (control)
% It is a matrix [3x3n x 1]
% N: number of nodes points
%
%OUTPUT
% FunVal: value of the function to minimize, in this case sum of the square
% of the control vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the function of which to find the root

% Extract state and control from y
u = y(7:9,1:N);

% Calculate the value of the function to minimize with CS quadrature
tVec = scalingFunction(tVec,1/(tVec(end)-tVec(1)));
tStep = tVec(3) - tVec(1);
uSquared = u(1,:).^2 + u(2,:).^2 + u(3,:).^2;
DV = (tStep/6)*(uSquared(1) + 2*sum( uSquared(3:2:end-1) ) + ...
         + 4*sum(uSquared(2:2:end-1)) + uSquared(end));
FunVal = 0.5*DV;

% Calculate the value of the function to minimize with just a summation
% uSquared = zeros(N,1);
% for i=1:N
%   uSquared(i) = norm(u(1:3,i))^2;
%   %uSquared(i) = norm(u(1:3,i));
% end
% FunVal = 0.5*sum(uSquared);

end

