function Kep = Pseudo2Kep(Pseudo)
% 
% Pseudo2Kep: Function to get Keplerian elements from Pseudo elements
%
% INPUT
% Pseudo: matrix with the Pseudo elements

%
% OUTPUT: 
% Kep: matrix with the keplerian elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conversion

% Should be a column vector
if size(Pseudo,1)==6
elseif size(Pseudo,2)==6
    Pseudo=Pseudo';
else
    error('Invalid input dimension')
end

% Extract Pseudo elem
p=Pseudo(1,:);
f=Pseudo(2,:);
g=Pseudo(3,:);
h=Pseudo(4,:);
k=Pseudo(5,:);
L=Pseudo(6,:);

% Obtain Kep elem
e=sqrt(f.^2+g.^2);
a=p./(1-e.^2);
i=2*atan(sqrt(h.^2+k.^2));

if i==0
    Om=0;
else
    Om=atan2(k./tan(i/2),h./tan(i/2));
end

if e==0
    om_s=0;
else
    om_s=atan2(g./e,f./e);
end

om=om_s-Om;
while Om<0
    Om=Om+2*pi;
end
while Om>=2*pi
    Om=Om-2*pi;
end

while om<0
    om=om+2*pi;
end
while om>=2*pi
    om=om-2*pi;
end

% M=l-om_s;
th=L-om_s;
while th<0
    th=th+2*pi;
end
while th>=2*pi
    th=th-2*pi;
end

% Gather elements in a matrix
Kep=[a; e; i; Om; om; th];

return