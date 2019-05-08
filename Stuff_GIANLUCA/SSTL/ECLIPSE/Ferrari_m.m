function Roots=Ferrari_m(Coeff)

alfa=-3*Coeff(2)^2/(8*Coeff(1)^2)+Coeff(3)/Coeff(1);
beta=Coeff(2)^3/(8*Coeff(1)^3)-Coeff(2)*Coeff(3)/(2*Coeff(1)^2)+Coeff(4)/Coeff(1);
gamma=-3*Coeff(2)^4/(256*Coeff(1)^4)+Coeff(2)^2*Coeff(3)/(16*Coeff(1)^3)-Coeff(2)*Coeff(4)/(4*Coeff(1)^2)+Coeff(5)/Coeff(1);
P=-alfa^2/12-gamma;
Q=-alfa^3/108+alfa*gamma/3-beta^2/8;
R=-Q/2+sqrt(Q^2/4+P^3/27);
U=(R)^(1/3);
if U==0
    y=-5/6*alfa-(Q)^(1/3);
else
    y=-5/6*alfa+U-P/(3*U);
end
W=sqrt(alfa+2*y);
D=sqrt(-(3*alfa+2*y+2*beta/W));
E=sqrt(-(3*alfa+2*y-2*beta/W));
Roots=[-0.25*Coeff(2)/Coeff(1)+0.5*(W+D);...
       -0.25*Coeff(2)/Coeff(1)+0.5*(W-D);...
       -0.25*Coeff(2)/Coeff(1)+0.5*(-W+E);...
       -0.25*Coeff(2)/Coeff(1)+0.5*(-W-E);];
   
return