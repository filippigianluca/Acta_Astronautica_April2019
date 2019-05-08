% Cristian Greco 4520319
%
% Function to convert from Keplerian parameters to Cartesian coordinates
%
% Inputs [m]/[-]/[rad]
% 
% Outputs: 
% r=3x1 [m]
% v=3x1 [m/s]

function [in]=Kep2Cart(Kep,mu)

if size(Kep,1)==6
    n = size(Kep,2);
elseif size(Kep,2)==6
    Kep = Kep';
    n = size(Kep,2);
else
    error('Invalid input dimension')
end
a=Kep(1,:);
e=Kep(2,:);
i=Kep(3,:);
OM=Kep(4,:);
om=Kep(5,:);
teta=Kep(6,:);

% Check on circular or elliptical orbit
if e>0 % Vectorial check
    if e<1
        
        % Semi-latus rectum
        p=a.*(1-e.^2);
        
        % Angular momentum
        H=sqrt(mu.*p);
        
        % Conversions
        l1 = cos(OM).*cos(om) - sin(OM).*sin(om).*cos(i);
        l2 = -cos(OM).*sin(om) - sin(OM).*cos(om).*cos(i);
        m1 = sin(OM).*cos(om) + cos(OM).*sin(om).*cos(i);
        m2 = -sin(OM).*sin(om) + cos(OM).*cos(om).*cos(i);
        n1 = sin(om).*sin(i);
        n2 = cos(om).*sin(i);
        
        % Perifocal coordinates
        r_pf = p./(1 + e.*cos(teta));
        x_pf = r_pf.*cos(teta);
        y_pf = r_pf.*sin(teta);
        z_pf = zeros(1,n);
        vett_r_pf = [x_pf ; y_pf ; z_pf];
        r = [];
        v = [];
        
        for i=1:n
            A = [ l1(i) l2(i) 0
                m1(i) m2(i) 0
                n1(i) n2(i) 0];
            
            r = [r,A*vett_r_pf(:,i)];
            
            
            % Velocities
            B = [-l1(i)*sin(teta(i)) + l2(i)*(e(i)+cos(teta(i)));
                -m1(i)*sin(teta(i)) + m2(i)*(e(i)+cos(teta(i)));
                -n1(i)*sin(teta(i)) + n2(i)*(e(i)+cos(teta(i)));];
            
            v = [v,(mu/H(i))*B];
        end
    else
        error('Not elliptical orbit');
    end
else
    error('Not elliptical orbit');
end
in = [r;v];
end