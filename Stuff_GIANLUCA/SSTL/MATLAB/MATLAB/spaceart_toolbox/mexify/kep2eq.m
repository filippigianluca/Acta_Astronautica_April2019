function Equin=kep2eq(Kep)

if size(Kep,1)==6
elseif size(Kep,2)==6
    Kep=Kep';
else
    error('Invalid input dimension')
end

a=Kep(1,:);
e=Kep(2,:);
i=Kep(3,:);
Om=Kep(4,:);
om=Kep(5,:);
th=Kep(6,:);

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
P1=e.*sin(om_s);
P2=e.*cos(om_s);
Q1=tan(i/2).*sin(Om);
Q2=tan(i/2).*cos(Om);

E=2*atan(sqrt((1-e)./(1+e)).*tan(th/2));
M=E-e.*sin(E);
% l=om_s+M;
% while l>=2*pi
%     l=l-2*pi;
% end
% while l<0
%     l=l+2*pi;
% end
L=th+om_s;
while L>=2*pi
    L=L-2*pi;
end
while L<0
    L=L+2*pi;
end

Equin=[a; P1; P2; Q1; Q2; L];

return