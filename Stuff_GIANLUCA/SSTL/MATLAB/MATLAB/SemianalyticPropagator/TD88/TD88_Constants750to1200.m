function [H,c1,p,knj,kn]=TD88_Constants750to1200()
% Values of constants for 750km to 1200km

% height scale (m)
H=27852.29; % meters

% numerical constants, "an" in paper
c1=[0.007243 0.1778 0.1449 -0.01179 7,011 6.968 3.301 14.91];

% phases, "pn" in paper
p=[0, 0, -263, 263, -29.41, 8.0913, 10.0813]*pi/180;

% numerical constants
knj=[
0.766348e-8 0.15645e-9 0.1943e-10
-0.440146e-8 0.36649e-9 0.2661e-9
0.118107e-9 -0.14007e-9 -0.1470e-11
-0.159664e-10 -0.65029e-11 -0.2317e-11
-0.240756e-9 -0.14318e-10 -0.3506e-11
0.643785e-10 0.14922e-9 0.2019e-10
0.744666e-11 0.44938e-11 0.2066e-11]; 

% numerical constants 
kn=[0.352814e-14 0.807687e-14 0.845184e-15 -0.116620e-16 -0.260648e-15 0.991392e-15 0.699853e-16];

end