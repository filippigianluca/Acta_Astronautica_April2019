% plot Belief 3Dload
clear all; close all; clc

load('decomposition_IAC_3_SEPTEMBER_mass')
F_mass = LIST.F_Bel;
Bel_mass = LIST.Bel;

load('decomposition_IAC_3_SEPTEMBER_data_volume')
F_DV = LIST.F_Bel;
Bel_DV = LIST.Bel;


for i=1:length(F_mass)
    for j=1:length(F_DV)
        
%         F3D(i,j) = [F_mass(i), F_DV(j)];
        BELIEF3D(i,j) = prod([Bel_mass(i) Bel_DV(j)]);
    end
end


[X,Y] = meshgrid(F_mass,F_DV);
surf(X, Y, BELIEF3D)
