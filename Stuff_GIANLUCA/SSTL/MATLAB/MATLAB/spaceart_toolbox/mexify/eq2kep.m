function Kep=eq2kep(Equin)

if size(Equin,1)==6
elseif size(Equin,2)==6
    Equin=Equin';
else
    error('Invalid input dimension')
end

a=Equin(1,:);
P1=Equin(2,:);
P2=Equin(3,:);
Q1=Equin(4,:);
Q2=Equin(5,:);
L=Equin(6,:);

e=sqrt(P1.^2+P2.^2);
i=2*atan(sqrt(Q1.^2+Q2.^2));
if i==0
    Om=0;
else
    Om=atan2(Q1./tan(i/2),Q2./tan(i/2));
end
if e==0
    om_s=0;
else
    om_s=atan2(P1./e,P2./e);
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

Kep=[a; e; i; Om; om; th];

return