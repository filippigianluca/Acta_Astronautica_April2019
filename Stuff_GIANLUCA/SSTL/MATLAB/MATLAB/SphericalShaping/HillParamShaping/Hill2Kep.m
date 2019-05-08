function Kep = Hill2Kep(Hill,mu)
% 
% Hill2Kep: Function to get Keplerian elements from Hill elements
%
% INPUT
% Hill: matrix with the Pseudo elements

%
% OUTPUT: 
% Kep: matrix with the keplerian elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion

% Should be a column vector
if size(Hill,1)==6
elseif size(Hill,2)==6
    Hill=Hill';
else
    error('Invalid input dimension')
end

% Extract Pseudo elem
r=Hill(1,:);
p=Hill(2,:);
G=Hill(3,:);
u=Hill(4,:);
H=Hill(5,:);
h=Hill(6,:);

% Obtain Kep elem
a = (G.^2.*r.^2*mu)./((mu*r).^2 - (G.^2 - mu*r).^2 - (G.*p.*r).^2);

e = sqrt( ((G.^2 - mu*r)./(mu.*r)).^2 + (G.*p/mu).^2);

sinTh = p.*(a.*(1-e.^2)/mu)./e;
cosTh = ((a.*(1-e.^2))./r - 1)./e;

th = atan2(sinTh,cosTh);

while th<0
    th=th+2*pi;
end
while th>=2*pi
    th=th-2*pi;
end

i=acos(H./G);
if i==0
    Om=0;
else
    Om=h;
end

while Om<0
    Om=Om+2*pi;
end
while Om>=2*pi
    Om=Om-2*pi;
end

om = u - th;

while om<0
    om=om+2*pi;
end
while om>=2*pi
    om=om-2*pi;
end

% Gather elements in a matrix
Kep=[a; e; i; Om; om; th];

return