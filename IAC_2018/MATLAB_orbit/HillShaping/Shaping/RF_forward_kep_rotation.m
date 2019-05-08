function kep_rot = RF_forward_kep_rotation(kep,omega0,angle)
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

A = acos2(-cos(B).*cos(C)+sin(B).*sin(C).*cos(a),sin(a)) ;

b = acos2((cos(B)+cos(C).*cos(A))./(sin(C).*sin(A)),sin(B));
c = acos2((cos(C)+cos(B).*cos(A))./(sin(B).*sin(A)),sin(C));

kep_rot(1:2,:) = kep(1:2,:) ;
kep_rot(3,:) = wrapTo2Pi(pi-A) ; 
kep_rot(4,:) = wrapTo2Pi(b+omega0) ;
kep_rot(5,:) = wrapTo2Pi(c + kep(5,:)) ;
kep_rot(6,:) = kep(6,:);

end