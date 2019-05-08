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
% Number of rows should be 6

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

%--------------------------------------------------------------------------
% Obtain Hill elements
%--------------------------------------------------------------------------
% r (radius)
pp = a.*(1 - e.^2);
r = pp./(1 + e.*cos(th));

% v_r (satellite range rate = projection of velocity in radius direction)
v_r = sqrt(mu./pp).*e.*sin(th);

% G (angular momentum)
G = sqrt(mu.*pp);

% u (anomaly from descending - or ascending ?? - node)
u = wrapTo2Pi( th + om );

% H (angular momentum normal to eq plane)
% H = G.*cos(i);

% h (longitude of ascending node)
h = wrapTo2Pi(Om); 

%--------------------------------------------------------------------------
% Gather Hill elements in a column vector (or matrix)
%--------------------------------------------------------------------------
Hill = [r; v_r; G; u; i; h];

return