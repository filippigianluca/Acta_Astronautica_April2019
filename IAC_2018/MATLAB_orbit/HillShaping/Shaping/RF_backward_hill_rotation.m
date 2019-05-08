function Hill = RF_backward_hill_rotation(Hill_rot,omega0,angle)
% Hold for angle < pi/2 I guess
% Use spherical trigonometry
if size(Hill_rot,2) == 6
    Hill_rot = Hill_rot';
elseif size(Hill_rot,1) == 6
else
    error('Wrong dimension')
end

c = Hill_rot(6,:) - omega0 + 1e-6 ;
if min(mod(c,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
    keyboard
end

A = angle ; 
if min(mod(A,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
    keyboard
end

B = pi - Hill_rot(5,:) ;
if min(mod(B,pi)) == 0
    error('Singular case for angle-side-angle spherical triangle')
    keyboard
end

C = acos2(-cos(A).*cos(B)+sin(A).*sin(B).*cos(c),sin(c)) ;

a = acos2((cos(A)+cos(B).*cos(C))./(sin(B).*sin(C)),sin(A));
b = acos2((cos(B)+cos(A).*cos(C))./(sin(A).*sin(C)),sin(B));

Hill(1:3,:) = Hill_rot(1:3,:) ;
Hill(5,:) = wrapTo2Pi(C) ; 
Hill(6,:) = wrapTo2Pi(b + omega0) ;
Hill(4,:) = wrapTo2Pi( Hill_rot(4,:) - a) ;

end