function [x]=eq2cart_m(Equin,mu)

% (c)  Massimiliano Vasile
% Modified by F.Zuiani (2010): 
% -Calculation of x(7:12) eliminated.
% -Inverted position of f,g and h,k in the input vector to comply with a
% different standard.
% -eq(1) is now the semi-major axis.
a=Equin(1);
f=Equin(3);
g=Equin(2);
h=Equin(5);
k=Equin(4);
L=Equin(6);
p=a*(1-f^2-g^2);

csi=1+f*cos(L)+g*sin(L);
rm=p/csi;
s2=1+h^2+k^2;
alf2=h^2-k^2;

r(1)=(rm/s2)*(cos(L)+alf2*cos(L)+2*h*k*sin(L));
r(2)=(rm/s2)*(sin(L)-alf2*sin(L)+2*h*k*cos(L));
r(3)=(2*rm/s2)*(h*sin(L)-k*cos(L));

v(1)=-(1/s2)*sqrt(mu/p)*(sin(L)+alf2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+alf2*g);
v(2)=-(1/s2)*sqrt(mu/p)*(-cos(L)+alf2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+alf2*f);
v(3)=(2/s2)*sqrt(mu/p)*(h*cos(L)+k*sin(L)+f*h+g*k);

x(1:3,1)=r;
x(4:6,1)=v;

return