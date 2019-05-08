function kep_rot = RF_forward_rotation(kep,omega0,angle)
% Hold for angle < pi/2 I guess
% Use spherical trigonometry
if size(kep,2) == 6
    kep = kep';
elseif size(kep,1) == 6
else
    error('Wrong dimension')
end

a = kep(4,:)-omega0 ; 
if min(mod(a,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end
C = angle ; 
if min(mod(C,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end
B = kep(3,:) ;
if min(mod(B,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end

cotb = (cos(a).*cos(C)+cot(B).*sin(C)) ./ sin(a) ;
b = acot(cotb) ;

cotc = (cos(a).*cos(B)+cot(C).*sin(B)) ./ sin(a) ;
c = acot(cotc) ;



sinA = sin(a).*sin(B)./sin(b) ;
cosA = (cos(a)-cos(b).*cos(c))./(sin(b).*sin(c));

A = atan2(sinA,cosA);

kep_rot(1:2,:) = kep(1:2,:) ;
kep_rot(3,:) = wrapTo2Pi(pi-A) ; 
kep_rot(4,:) = wrapTo2Pi(b+omega0) ;
kep_rot(5,:) = wrapTo2Pi(c + kep(5,:)) ;
kep_rot(6,:) = kep(6,:);

end