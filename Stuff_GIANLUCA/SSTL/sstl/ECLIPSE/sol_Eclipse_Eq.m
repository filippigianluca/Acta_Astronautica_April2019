function [L_i, L_o, kep, dt, kep2, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt, mu, R)

% Inputs
% Eq: current set of Equinoctial elements in the J2000, Earth-centered equatorial frame.
% tt: current time in MJD2000 (in days).
% mu: gravitational constant of the Earth.
% R: planetary radius of the Earth for eclipse calculation.

% Outputs
% L_i: True longitude of the eclipse entry point.
% L_o: True longitude of the eclipse exit point.
% kep: keplerian parameters in the J2000, Earth-centered equatorial frame,
% corresponding to Eq
% dt: eclipse duration.
% kep2: keplerian parameters in an Earth-centered, Ecliptic reference frame
% in which the x direction is aligned with the Earth-Sun vector.
% alpha_Sun: azimuth of the Sun-Earth direction (i.e. of the solar
% radiation pressure acceleration vector) in the r-t-h reference frame.
% beta_Sun: elevation of the Sun-Earth direction (i.e. of the solar
% radiation pressure acceleration vector) in the r-t-h reference frame.
% phi_Sun:  ?
%
% Federico Zuiani 2012-2013
% Comments and code tyding: Marilena Di Carlo, 2015 (comments may be
% wrong!)

% Note: all units in km, s, etc. Or at least properly rescaled.
 
% tt is the time in MJD2000!


% Earth position vector at time tt in the heliocentric reference frame [km]
rr = EphSS(3,tt);

% Earth azimuth angle at time tt [rad]
az_S_E = atan2(rr(2), rr(1));

% Keplerian elements of the point given in equinoctial elements as input
% (unit are the same as Eq in input) - in the Earth centered reference
% frame
kep = eq2kep(Eq);

% If also the angles are required as output
if nargout > 5

    % Ecliptic inclination [rad]
    incl_ecl = astro_constants(8);

    % Rotate rr: from ecliptic to equatorial reference place
    % Clockwise rotation of an angle equal to the ecliptic
    % inclination around the x axis
    % The rotation matrix would be 
    % R = [1 0 0; 0 cosd(-23) sind(-23); 0 -sind(-23) cosd(-23)]
    % The terms multiplying rr(3) have been negletected because rr(3)=0
    % (the Earth IS on the plane of the ecliptic)
    r_S_E_eq = [rr(1); cos(incl_ecl)*rr(2); sin(incl_ecl)*rr(2)];
% keyboard
    % Transformation to radial-transversal-h reference frame
    r_S_E_rth = car_rthT(r_S_E_eq, kep2cart(kep,mu));

    rr = norm(r_S_E_rth);
    
    % Sun elevation angle wrt spacecraft
    beta_Sun=asin(r_S_E_rth(3)/rr);
    
    % Sun azimuth angle wrt spacecraft
    alpha_Sun=atan2(r_S_E_rth(2)/rr/cos(beta_Sun),r_S_E_rth(1)/rr/cos(beta_Sun));

    % If also phi_Sun is required as output:
    if nargout > 7
        
        % Transfrom vector from Sun to Earth expressed in Earth reference
        % frame into rth reference frame of satellite at perigee
        r_S_E_eth = car_rthT(r_S_E_eq, kep2cart([kep(1:5),0],mu));
        
        % Angle computed in counter clockwise direction from Sun-Earth
        % vector in rth reference frame to r direction of this reference
        % frame
        phi_Sun = mod(2*pi-atan2(r_S_E_eth(2),r_S_E_eth(1)),2*pi);
    end
end

% Compute polynomial coefficients
[Coeff, kep2] = Eclipse_Eq(kep, az_S_E, mu, R);

if any(isnan(Coeff)|isinf(Coeff))||kep(2)>=1
    warning('Some element of Coeff is Nan or Inf!')
    format longg
    Eq
    tt
    mu
    R
    L_i=NaN;
    L_o=NaN;
    kep=NaN;
    kep2=NaN;
    dt=NaN;
    return
end


% Find polynomial roots
%% Numerical

% sol = roots(Coeff);

%% Ferrari's Method

sol = Ferrari_m(Coeff);

%% Identify correct solutions

th = [];
c_th = zeros(1,length(sol));

% -------------------------------------------------------------------------
% Commented by Federico:
% th_x=NaN*ones(2,length(sol));
% -------------------------------------------------------------------------

% I guess that some of the solution may have imaginary component and here
% Federico solves the problem... don't know how...
for i = 1 : length(sol)
    
    if ( abs(real(sol(i))) < 1 + 1e-6)   &&   ( (abs(real(sol(i))) > abs(imag(sol(i)))) || (abs(imag(sol(i)))<=1e-7) )
        p_flag = 0;
        c_th(i) = real(sol(i));
        % +sin
        s_th = sqrt(1-c_th(i)^2);
        if any(~isreal([s_th c_th(i)]))
            s_th=real(s_th);
        end
        
        kep2(6)=atan2(s_th,c_th(i));
        x=kep2cart(kep2,mu);
        c_1_p=x(1);
        c_2_p=sqrt(x(2)^2+x(3)^2)-R;
%         fprintf('%d\t%d\t%d\n',kep2(6),c_1_p/R,c_2_p/R)
        if (c_1_p/R>-1e-3)&&(abs(c_2_p)/R<1e-3)
            th=[th kep2(6)];
            p_flag=1;
%             th_x(1,i)=kep2(6);
        end
        % -sin
        s_th=-s_th;

        kep2(6)=2*pi-kep2(6);
        x=kep2cart(kep2,mu);
        c_1_m=x(1);
        c_2_m=sqrt(x(2)^2+x(3)^2)-R;
%         fprintf('%d\t%d\t%d\n',kep2(6),c_1_m/R,c_2_m/R)
        if (c_1_m/R>-1e-3)&&(abs(c_2_m)/R<1e-3)
            if p_flag
                if ((abs(c_2_m)<abs(c_2_p)/100)||(c_1_m>100*c_1_p))
                    th(end)=kep2(6);
                end
            else
                th=[th kep2(6)];
            end
%             th_x(2,i)=kep2(6);
        end

    end
    
end

% If the previous correction to identify correct solution did return an
% empty th vector, then there are no eclipses
if isempty(th)
    % -------------------------------------------------------------------------
    % Commented by Federico:
        % 
    %     %%
    %     nn=100;
    %     c_th_t=linspace(-1,1,nn);
    %     ff=zeros(1,nn);
    %     for i=1:nn
    %         ff(i)=f_ecl(Coeff,c_th_t(i));
    %     end
    %     hh=figure(23);
    %     plot(c_th_t,ff,'b-',cos(th),0*th,'r*',c_th,0*c_th,'go')
    %     axis([-1 1 -Inf Inf])
    %     grid on
    %     keyboard
    %     close(hh)
    %     %%
    % -------------------------------------------------------------------------
    L_i  = NaN;
    L_o  = NaN;
    kep  = NaN;
    kep2 = NaN;
    dt   = NaN;
    return
end

th = sort(th);
th = [th(1) th(end)];
th = mod(th,2*pi);

% Identify eclipse entrance and exit
while th(2)-th(1)<0
    th(2)=th(2)+2*pi;
end

while th(2)-th(1)>=2*pi
    th(2)=th(2)-2*pi;
end

th_m=(th(1)+th(2))/2;
kep2(6)=th_m;
x=kep2cart(kep2,mu);
c_1=x(1);
c_2=sqrt(x(2)^2+x(3)^2)-R;

if (mod(th(2)-th(1),2*pi)<=pi)&&(c_1/R>-1e-3)&&(c_2/R<1e-3)
    L_i=th(1)+kep(4)+kep(5);
    L_o=th(2)+kep(4)+kep(5);
else
    L_i=th(2)+kep(4)+kep(5);
    L_o=th(1)+kep(4)+kep(5)+2*pi;
    th=[th(2) th(1)];
end

while L_o-Eq(6)<=0
    L_i=L_i+2*pi;
    L_o=L_o+2*pi;
end

while L_o-Eq(6)>2*pi
    L_i=L_i-2*pi;
    L_o=L_o-2*pi;
end

if th(2)-th(1)>pi
    th=[th(2) th(1)];
    save('Ze_error','Eq','tt','mu','R')
    warning('Eclipse arc spans more than 180 degrees!')
end
if nargout>3
    if th(2)<th(1)
        th(2)=th(2)+2*pi;
    end
        
    % Find the time corresponding to a certain true anomaly
    % INPUT:
    %   1st:    True anomaly [rad].
    %   2nd:    Semi-major axis [L].
    %   3rd:    Eccentricity.
    %   4th:   Planetary constant (mu = mass * G) [L^3/T^2].
    %   5th:   Fixed true anomaly [rad].
    %   6th:   Time corresponding to f0 [T].
    %           > Optional: Default value = 0.
    
    
    dt = kepEq_t_matlab(th(2), kep2(1), kep2(2), mu, th(1), 0);
%         dt = kepEq_t(th(2), kep2(1), kep2(2), mu, th(1), 0);

    kep2(6)=kep(6);
end


% fprintf('%d\t%d\n',L_i,L_o)

return

% -------------------------------------------------------------------------
% Commented by Federico:
% function f=f_ecl(Coeff,x)
% 
% f=sum(Coeff.*x.^[(length(Coeff)-1):-1:0]');
% 
% return
% -------------------------------------------------------------------------