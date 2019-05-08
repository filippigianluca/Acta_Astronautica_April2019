function Pseudo = Kep2Pseudo(Kep)
% 
% Kep2Pseudo: Function to get Pseudo elements from Keplerian elements
%
% INPUT
% Kep: matrix with the keplerian elements
% mu: gravitational parameter
%
% OUTPUT: 
% Pseudo: matrix with the Pseudo elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion

% Should be a column vector
if size(Kep,1)==6
elseif size(Kep,2)==6
    Kep=Kep';
else
    error('Invalid input dimension')
end

% Extract Kep elem
a=Kep(1,:);
e=Kep(2,:);
i=Kep(3,:);
Om=Kep(4,:);
om=Kep(5,:);
th=Kep(6,:);

% Obtain Pseudo elem
om_s=Om+om;
while om_s>=2*pi
    om_s=om_s-2*pi;
end
while om_s<0
    om_s=om_s+2*pi;
end
while Om>=2*pi
    Om=Om-2*pi;
end
while Om<0
    Om=Om+2*pi;
end

p = a.*(1 - e.^2);
f=e.*cos(om_s);
g=e.*sin(om_s);
h=tan(i/2).*cos(Om);
k=tan(i/2).*sin(Om);

L=th+om_s;
while L>=2*pi
    L=L-2*pi;
end
while L<0
    L=L+2*pi;
end

% Gather elements in a matrix
Pseudo=[p; f; g; h; k; L];

return