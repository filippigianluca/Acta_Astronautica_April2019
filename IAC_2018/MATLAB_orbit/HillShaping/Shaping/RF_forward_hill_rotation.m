function Hill_rot = RF_forward_hill_rotation(Hill,omega0,angle)
% Hold for angle < pi/2 I guess
% Use spherical trigonometry
if size(Hill,2) == 6
    Hill = Hill';
elseif size(Hill,1) == 6
else
    error('Wrong dimension')
end

a = Hill(6,:) - omega0 ; 
if min(mod(a,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end
C = angle ; 
if min(mod(C,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end
B = Hill(5,:) ;
if min(mod(B,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
end

A = acos2(-cos(B).*cos(C)+sin(B).*sin(C).*cos(a),sin(a)) ;

b = acos2((cos(B)+cos(C).*cos(A))./(sin(C).*sin(A)),sin(B));
c = acos2((cos(C)+cos(B).*cos(A))./(sin(B).*sin(A)),sin(C));

Hill_rot(1:3,:) = Hill(1:3,:) ;
Hill_rot(5,:) = wrapTo2Pi(pi-A) ; 
Hill_rot(6,:) = wrapTo2Pi(b+omega0) ;
Hill_rot(4,:) = wrapTo2Pi(c + Hill(4,:)) ;

end