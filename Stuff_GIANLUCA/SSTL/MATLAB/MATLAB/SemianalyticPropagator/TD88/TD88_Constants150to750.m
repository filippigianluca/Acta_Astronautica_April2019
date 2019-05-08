function [H,c1,p,knj,kn]=TD88_Constants150to750()
% Values of constants for 150km to 750km

% height scale (m)
H=29000; % meters

% numerical constants, "an" in paper
c1=[0.007 0.2875 0.04762 0.0471 7 7 0.3333 15];

% phases, "pn" in paper
p=[0, 0, -263, 263, -29.41, 8.0913, 10.0813]*180/pi;

% numerical constants
knj=[
0.766373e-8 0.16573e-9 0.3871e-10
-0.440149e-8 0.33428e-9 0.9352e-10
0.118107e-9 -0.14781e-9 -0.1518e-11
-0.159664e-10 -0.64670e-11 -0.2050e-11
-0.240755e-9 -0.13985e-10 -0.3059e-11
0.643785e-10 0.13618e-9 0.3517e-10
0.744666e-11 0.45416e-11 0.2080e-11];

% numerical constants 
kn=[0.296815e-14 0.281456e-13 -0.123300e-13 -0.114892e-16 -0.390065e-15 0.742439e-14 -0.341594e-15]; 
end
