function teta = M2teta(M,e)

% M2teta: Converts the mean anomaly into the true anomaly
%
% INPUTs: - M: Mean anomaly
%         - e: Eccentricity
% 
% OUTPUT: - teta: True anomaly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cristian Greco, 08-09-2016
% Formula taken by: Fundamentals of Astrodynamics, Wakker, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if e<1e-12
    om=0;
    teta=M; % teta=E + E=M --> teta=M
    display('Circular orbit')
else
    % Find E from M, iterative process
    tol=100;
    E=M; %Initial guess
    
    while tol>1e-15
        E_old=E;
        E=E+(M-E+e*sin(E))/(1-e*cos(E));
        tol=abs(E-E_old);
    end
    
    % Find teta from E, no ambiguity, E/2 and teta/2 alomays in the same
    % quadrant
    teta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    teta = wrapTo2Pi(teta);
end