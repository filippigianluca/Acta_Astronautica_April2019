function Hill = Kep2Hill(Kep,mu)
% 
% Kep2Hill: Function to get Hill elements from Keplerian elements
%
% INPUT
% Kep: matrix with the keplerian elements
% mu: gravitational parameter
%
% OUTPUT: 
% Hill: matrix with the Hill elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
% Number of column should be 6
if size(Kep,1)==6
elseif size(Kep,2)==6
    Kep = Kep';
else
    error('Invalid input dimension')
end

% Extract Kepler elem
a=Kep(1,:);
e=Kep(2,:);
i=Kep(3,:);
Om=Kep(4,:);
om=Kep(5,:);
th=Kep(6,:);

% Transform Om in the range [0-2pi]
while Om>=2*pi
    Om=Om-2*pi;
end
while Om<0
    Om=Om+2*pi;
end

%--------------------------------------------------------------------------
% Obtain Hill elements
%--------------------------------------------------------------------------
% r (radius)
pp = a.*(1 - e.^2);
r = pp./(1 + e.*cos(th));

% Obtain p
%[R,V] = KeplElem2rv(a,e,i,om,Om,th,mu);
%unitR = R/norm(R);
%p = dot(unitR,V);

% p (satellite range rate = projection of velocity in radius direction)
p = sqrt(mu./pp).*e.*sin(th);

% G (angular momentum)
G = sqrt(mu*pp);

% u (anomaly from descending - or ascending ?? - node)
u= th + om;
while u>=2*pi
    u=u-2*pi;
end
while u<0
    u=u+2*pi;
end

% H (angular momentum normal to eq plane)
H = G*cos(i);

% h (longitude of ascending node)
h = Om;
while h>=2*pi
    h=h-2*pi;
end
while h<0
    h=h+2*pi;
end

%--------------------------------------------------------------------------
% Gather Hill elements in a column vector (or matrix)
%--------------------------------------------------------------------------
Hill = [r; p; G; u; H; h];
return