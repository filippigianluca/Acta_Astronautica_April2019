function [dv,vin,error,xp2,vp2]=mimga3_2(t0,T,s,multirev,v0,options)
% mimga3 - Multi gravity assists, with multi revolution Lambert
%
%   This function computes the dv for a sequence of gravity assists and
%   deep space manoeuvres, in order to reach a given senquence of planets.
%   It uses the traditional patched-conic approximation, and the main body
%   during deep space flights is the Sun. It uses a multi revolution
%   lambert solver. Only flybys of planets (1 <= s <= 9) are allowed.
%
%   [dv, vout, error, vbrake, xp2, vp2] = mimga3_2(t0, T, s, v0, options)
%
% INPUT
%   t0      Departure time [d, MJD2000].
%   T       Vector of parameters.
%           If options(1)=0:
%               T = [alpha1, TOF1,
%                    gamma2, rp2, alpha2, TOF2,
%                    gamma3, rp3, alpha3, TOF3,
%                    ...];
%           If options(1)=1:
%               T = [gamma1, rp1, alpha1, TOF1,
%                    gamma2, rp2, alpha2, TOF2,
%                    ...];
%           If options(2)=0:
%               T = [...,
%                    gamman, rpn, alphan, TOFn];
%           If options(2)=1:
%               T = [...,
%                    gamman, rpn, alphan, TOFn,
%                    gamma, rp];
%           Where:
%               TOF: Time of flight [d];
%               gamma: angle of hyperbola plane [rad];
%               rp: radius of the pericentre of the hyperbola
%                   (adimensionalised with the radius of the planet).
%               alpha: fraction of TOF before the deep space manoeuvre.
%   v0      Initial velocity vector [km/s]. It can be relative to the first
%           planet or absolute, depending on options(3).
%   s       Planetary sequence. The first planet is the departure planet,
%           the last planet is the arrival planet. The flyby of the
%           starting and the ending planet is regulated by options(1:2).All
%           other planets will be flown by. Depending of the options, s
%           must contain at least one or two planets.
%   multirev    Matrix containing, in each column (relative to each Lambert
%               arc), the type of transfer (direct or retrograde), the
%               number of revolutions for the Lambert arc, the case of
%               transfer (low or high energy):
%                   multirev = [type1, type2,  ...;
%                               nrev1, nrev2, ...;
%                               case1, case2, ...];
%               Leave it empty [] for using 0 rev, low energy, direct
%               transfer for each arc. Or it's possible to fill only the
%               first 1 or 2 rows, to take 0 rev and low energy.
%   options Vector of options. It is optional, and it is also possible to
%           specify only the first options on the left of the vector. All
%           the others are assumed 0.
%           options(1) determines how to start the flight.
%               0:  The spacecraft is at time t0 at the first
%                   planet with velocity v0, but no flyby is performed, so
%                   the path starts with a deep space flight, then a deep
%                   space manoeuvre, a flyby of the next planet and so on,
%                   until the last planet. This case can be used for the
%                   launch case, in which alpha1=0, so that the deep space
%                   manoeuvre is the launch.
%               1:  The spacecraft is at time t0 at the first
%                   planet with velocity v0. It performs a flyby, then a
%                   deep space manoeuvre, and so on, until the last planet.
%           options(2) determines how to end the flight at the last planet.
%               0:  The function ends at the encounter of the last planet,
%                   without doing the flyby of it.
%               1:  The function ends after a flyby of the last planet.
%           options(3) is used to know whether v0 is absolute or relative
%               to the first planet of the sequence.
%               0:  v0 is the absolute initial velocity of the spacecraft,
%                   given in cartesian coordinates.
%               1:  v0 input is relative to the first planet, in cartesian
%                   coordinates. In this case, if v0 = [0 0 0], it means
%                   that the spacecraft starts *on* the planet: so, if
%                   alpha1=0, is the launch case.
%               2:  v0 input is relative to the first planet, and given as
%                   the modulus, the in-plane angle with respect to the
%                   planetary velocity, and the out-of-plane angle [rad].
%               3:  As in case 2, but the 2 angles are given adimensionally
%                   in the interval [0, 1], such that the sphere is sampled
%                   uniformly. See Weisstein, Eric W. "Sphere Point
%                   Picking." From MathWorld--A Wolfram Web Resource.
%                   http://mathworld.wolfram.com/SpherePointPicking.html
%           options(4) is used to plot the trajectory or print data.
%               0: No plot.
%               1: The trajectory is plotted in a new figure.
%               2: Plot and print trajectory data on the screen.
%               
% OUTPUT
%   dv      Sequence of dv [km/s]. The first element is the launch dv (when
%           available, i.e. options(1)=0). The last element is the deep
%           space dv before encountering the plast planet. Elements between
%           them are deep space dvs. No brake dv is computed.
%   vout    Absolute velocity vector at the end of the sequence [km/s].
%   error   0 if and only if everything ok. If error occurs, dv = 1e6.
%   xp      Vector of the position of the last encountered planet [km].
%   vp      Vector of the velocity of the last encountered planet [km/s]
%
% EXAMPLES
%   Launch from Earth and deep space flight to Venus:
%       t0 = [3359.38387247647];
%       T = [0 170.49933001840];
%       s = [3,2];
%       [dv,vout,error,xp2,vp2] = mimga3_2(t0, T,s,[],[0,0,0],[0,0,1]);
%   
%   Just a flyby of Venus:
%       t0 = [3359.38387247647+170.49933001840];
%       T = [3.26576512201 2.07228045370];
%       s = [2];
%       vin = [... ... ...];
%       [dv,vout,error,xp2,vp2] = mimga3_2(t0, T,s,[],vin,[1,0,0]);
%
%   Launch from Earth, deep space flight and flyby of Venus:
%       t0 = [3359.38387247647];
%       T = [0 170.49933001840 3.26576512201 2.07228045370];
%       s = [3,2];
%       [dv,vout,error,xp2,vp2] = mimga3_2(t0, T,s,[],[0,0,0],[0,1,1]);
%
%   Complete EVEEJ transfer, with plot:
%       t0 = [3359.38387247647];
%       T = [0, 170.49933001840,...
%            3.26576512201, 2.07228045370, 0.06197382751, 320.21104182897, ...
%            2.65131587336, 1.66668311162, 0.07728172747, 730.48301546338, ...
%            3.16200529828, 1.41291334360, 0.10740729409, 847.09714061671];
%       s = [3,2,3,3,5];
%       [dv,vout,error,xp2,vp2] = mimga3_2(t0,T,s,[],[0,0,0],[0,0,1,1]);
%
% CALLED FUNCTIONS
%   EphSS, lambertMR, swingby, keppro3, astro_constants, d2b_eq,
%   cart2kep.
%
% Freely inspired by Massimiliano Vasile's mimga.
% Matteo Ceriotti, 11-11-2006
%
% Modified: Matteo Ceriotti, 01-08-2006: Added support for inital velocity
%           given with modulus and 2 angles.
% Modified: Matteo Ceriotti, 09-11-2006: Does not read the radius and the
%           mu of the planet anymore, if no flyby will be performed.
% Modified: Matteo Ceriotti, 10-11-2006: When plotting, chooses the color
%           as a function of the semimajor axis of the planet
% Version 3 - Matteo Ceriotti, 01-01-2007: Multi rev. Lambert
% Modified: Matteo Ceriotti, 13-02-2007: Uses keppro3, cart2kep.
% Modified: Matteo Ceriotti, 23-07-2007: When error occurs, dv = 1e6.
% Modified: Matteo Ceriotti, 26-07-2007: Added verbose possibility in
%           options(4).
% Modified: Matteo Ceriotti, 30-07-2007: Added the possibility of uniform
%           sampling of the sphere for v0.
% Modified: Matteo Ceriotti, 02-02-2009: Changed two_body_dynamics to
%           d2b_eq. Fixed row/column vector problem after swingby call.
%
% -------------------------------------------------------------------------

if nargin==5
    options=[];
end

if options(4)==2 % Verbose
    fprintf('  *** mimga3_2 - Verbose mode ***\n');
    fprintf('  ***\nInput parameters\n');
    fprintf('t0 = %f\nT =\n[',t0);
    for i=1:length(T)
        fprintf('%f ',T(i));
    end
    fprintf(']\ns =\n[ ');
    for i=1:length(s)
        fprintf('%d ',s(i));
    end
    fprintf(']\nmultirev = [\n');
    for j=1:size(multirev,1)
        for i=1:size(multirev,2)
            fprintf('%d ',multirev(j,i));
        end
        fprintf(';\n');
    end
    fprintf('v0 =\n[');
    for i=1:length(v0)
        fprintf('%f ',v0(i));
    end
    fprintf(']\noptions =\n[');
    for i=1:length(v0)
        fprintf('%d ',options(i));
    end
    fprintf(']\n');
end

% Completes options vector
if length(options)<4
    options=[options,zeros(1,4-length(options))];
end
if options(4) % Plot flag
    % Initialises figure
    hold on;
    line(0,0,0,'Marker','*');
    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    AU=astro_constants(2);
    % Initialise integrator used to plot trajectories
    ode45options = odeset('abstol',1e-9,'reltol',1e-7);
end
% Completes multirev matrix
nlambertarcs=length(s)-1;
if size(multirev,1)==0 % Not specified anything
    multirev(1,:)=zeros(1,nlambertarcs); % Assume direct (1) transfer for each arc
end
if size(multirev,1)==1 % Not specified the number of revolutions
    multirev(2,:)=zeros(1,nlambertarcs); % Assume 0 revolutions for each arc
end
if size(multirev,1)==2 % Not specified the case
    multirev(3,:)=zeros(1,nlambertarcs); % Assume low energy case (0) transfer for each arc
end

error=0;
mu=astro_constants(4); % Planetary constant of the Sun

planet1=s(1);
time=t0;
Tindex=1;
sindex=2;
[xp1,vp1]=EphSS(planet1,time); % Ephemeris of planet1

if options(4) % Plot flag
    line(xp1(1)/AU,xp1(2)/AU,xp1(3)/AU,'Marker','*') % Marker at departure point
end

dv=[];
switch options(3)
    case 0, % v0 is absolute, so nothing to do.

    case 1, % v0 is relative to the first planet, and cartesian
        v0=v0+vp1; % Conversion to absolute initial velocity
    case {2, 3}, % v0 is relative to the first planet, and given as modulus and 2 angles: the in-plane and the out-of plane.
        if options(3)==3 % A pre-transformation is needed, as the angles are given in adimensional components to uniformly sample the sphere
            % Uniform sampling of the sphere
            v0(2)=v0(2)*2*pi-pi/2;
            v0(3)=acos(2*v0(3)-1)-pi/2;
        end
        mod_v=v0(1); % Modulus of initial velocity relative to the departure planet
        psi=v0(2)-pi/2; % In-plane rotation
        theta=v0(3); % Out-of-plane rotation
        % Transposed Euler angles matrix in the special case phi=0.
        R_T=[cos(psi), -cos(theta)*sin(psi),  sin(theta)*sin(psi);...
             sin(psi),  cos(theta)*cos(psi), -sin(theta)*cos(psi);...
                0    ,       sin(theta)    ,       cos(theta)    ];
        v_rel_tnh=R_T*[0;mod_v;0]; % Relative velocity in the tangent-normal-binormal (h) reference frame
        v0=tnh_carT(v_rel_tnh,[xp1,vp1])';
        v0=v0+vp1; % Conversion to absolute initial velocity
    otherwise
        disp('Error: bad options(3)!');
        return
end

vp2=vp1; % Just as initialisation. If there is only one planet, this is useful for the calculation of the braking dv.
xp2=xp1; % Just as initialisation. If there is only one planet, this is useful for the output ephemeris.

switch options(1)
    case 0,
        % Same as (2), but the v0 velocity vector is after the flyby.
        % Basically, the flyby has been already done, and so starts the deep
        % space flight.
        vout=v0; % v0 is the velocity at the end of the flyby
    case 1,
        % Using the velocity vector v0 at the planet encounter, starts a swing
        % by of the first planet in the list, then continues with deep space
        % manoeuvre, swingby on the second planet, and so on (till the last
        % planet).
        vin=v0; % v0 is the velocity before the first flyby
    otherwise
        disp('Error: bad options(1)!');
        return
end

% Remeber to initialise Tindex, sindex, time, vin, planet1, xp1, vp1, nflybys, dv before the for

% Flybys and deep space manoeuvres
while Tindex<=length(T)
    % Here planet1 is the flyby planet, and planet2 is the planet towards
    % which the spacecraft is pointing (if available).
    
    % Flyby of planet1
    if sindex==2 && options(1)==0 % sindex==2 means that we are at the first planet
        % Skip the flyby of the first planet.
    else
        % Flyby of planet1
        % Reads mu and radius of the planet to be flown by.
        mup=astro_constants(10+planet1); % Gravity constant of the planet1
        rplanet=astro_constants(20+planet1); % Radius constant of the planet1

        % Reads flyby variables from vector T and s
        rp=T(Tindex+1)*rplanet; % Pericentre of the hyperbola (dimensional)
        psi=T(Tindex+0); % Hyperbola plane angle
        Tindex=Tindex+2;
        
        % Swing-by
        vinf1=vin-vp1; % Velocity relative to the planet
        n=cross(vinf1,vp1); % Reference vector
        n=n/norm(n);
        vinf2=swingby(vinf1,rp,mup,psi,n); % Swingby
        vout=vinf2(:)'+vp1(:)'; % Absolute velocity after the swingby
        if options(4)==2 % Verbose
            temp = car_rthT(vinf1,[xp1,vp1]);
            fprintf('  ***\nSwingby of %d - time = %f d,MJD2K\n',planet1,time);
            fprintf('Relative arrival velocity in r-th-h (also called R-S-T) [km/s]:\n[%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('|vinf| = %f km/s\n',norm(vinf1));
            fprintf('gamma = %f rad; rp = %f km (hp = %f km)\n',psi,rp,rp-rplanet);
            temp = car_rthT(vinf2,[xp1,vp1]);
            fprintf('Relative departure velocity in r-th-h (also called R-S-T) [km/s]:\n[%f %f %f]\n',temp(1),temp(2),temp(3));
        end
        
    end
    
    % Deep space flight (from planet1 to planet2) and manoeuvre
    if sindex>length(s) % sindex>length(s) means that we are at the last planet
        % Do nothing: only the flyby of the last planet had to be done.
        vin=vout;
    else
        % Reads trajectory variables from vector T and s
        alpha=T(Tindex+0); % Time fraction before deep space manoeuvre
        tof=T(Tindex+1); % Time of flight to the next planet
        planet2=s(sindex); % Number of second planet (arrival)

        % Propagation
        [sds,error]=keppro3([xp1,vout],tof*alpha*86400,mu); % Keplerian propagation
        xds=sds(1:3);vds1=sds(4:6);
        if error
            dv=1e6;
            vin=[];
            return
        end
        % Lambert arc
        [xp2,vp2]=EphSS(planet2,time+tof); % Position and velocity of the arrival planet
        [a_lambert,p_lambert,e_lambert,error,vds2,vin]=lambertMR(xds,xp2,tof*(1-alpha)*86400,mu,multirev(1,sindex-1),multirev(2,sindex-1),multirev(3,sindex-1),options(4));
        if error
            dv=1e6;
            vin=[];
            return
        end

        Tindex=Tindex+2; % Index of vector T
        sindex=sindex+1; % Next planet

        if options(4)==2 % Verbose
            fprintf('  ***\nDeep space flight %d-%d. time1 = %f d,MJD2K\n',planet1,planet2,time);
            fprintf('Initial absolute cartesian velocity [km/s]:\n[%f, %f, %f]\n',vout(1),vout(2),vout(3));
            temp = car_rthT(vout-vp1,[xp1,vp1]);
            fprintf('Initial relative velocity in r-th-h (also called R-S-T) [km/s]:\n[%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('Tof to the next planet = %f d; alpha to DSM = %f\n',tof,alpha);
            temp=cart2kep([xp1,vout],mu);
            fprintf('Orbital parameters before the DSM:\n');
            fprintf('a = %f km; e = %f; i = %f deg;\n',temp(1),temp(2),temp(3)/pi*180);
            fprintf('OM = %f deg; om = %f deg; th = %f deg;\n',temp(4)/pi*180,temp(5)/pi*180,temp(6)/pi*180);
            fprintf('rp = %f km; ra = %f km.;\n',temp(1)*(1-temp(2)),temp(1)*(1+temp(2)));
            fprintf('DSM time = %f d,MJD2K; dv = %f km/s\n',tof*alpha,norm(vds2-vds1));
            fprintf('DSM cartesian position [km]:\n[%f, %f, %f]\n',xds(1),xds(2),xds(3));
            temp=cart2kep([xds,vds2],mu);
            fprintf('Orbital parameters after the DSM:\n');
            fprintf('a = %f km; e = %f; i = %f deg;\n',temp(1),temp(2),temp(3)/pi*180);
            fprintf('OM = %f deg; om = %f deg; th = %f deg;\n',temp(4)/pi*180,temp(5)/pi*180,temp(6)/pi*180);
            fprintf('rp = %f km; ra = %f km.;\n',temp(1)*(1-temp(2)),temp(1)*(1+temp(2)));
            fprintf('Final absolute cartesian velocity [km/s]:\n[%f, %f, %f]\n',vin(1),vin(2),vin(3));
            temp = car_rthT(vin-vp2,[xp2,vp2]);
            fprintf('Final relative velocity in r-th-h (also called R-S-T) [km/s]:\n[%f %f %f]\n',temp(1),temp(2),temp(3));
            fprintf('Lambert parameters: D(0)/R(1) = %d; nrev = %d; Energy = %d\n',multirev(1,sindex-2),multirev(2,sindex-2),multirev(3,sindex-2));
        end

        % Plot
        if options(4) % Plot flag
            % Plots all the planets from time to time+tof
            for i=1:length(s)
                % Color of the line as a function of semimajor axis
                [x_color,v_color]=EphSS(s(i),0); % Ephemeris of s(i)
                temp=cart2kep([x_color,v_color],mu);
                a_color(i)=temp(1);
            end
            a_color=mod(a_color,(.1*AU))/(.1*AU);
            for i=1:length(s)
                [xpi,vpi]=EphSS(s(i),time); % Position and velocity of the arrival planet
                [ti,yi]=ode45(@d2b_eq,[0, tof*86400],[xpi,vpi],ode45options,mu);
                line(yi(:,1)/AU,yi(:,2)/AU,yi(:,3)/AU,'Color',[a_color(i), 1-a_color(i), 0]);
            end
            % Plots the spacecraft trajectory
            if alpha>0
                [tsp,ysp]=ode45(@d2b_eq,[0, tof*alpha*86400],[xp1,vout],ode45options,mu);
                line(ysp(:,1)/AU,ysp(:,2)/AU,ysp(:,3)/AU,'Color','k');
            end
            [tsp,ysp]=ode45(@d2b_eq,[0, tof*(1-alpha)*86400],[xds,vds2],ode45options,mu);
            line(ysp(:,1)/AU,ysp(:,2)/AU,ysp(:,3)/AU,'Color','b');
            line(ysp(1,1)/AU,ysp(1,2)/AU,ysp(1,3)/AU,'Marker','.') % Marker in dsm point
            line(xp2(1)/AU,xp2(2)/AU,xp2(3)/AU,'Marker','o') % Marker in rendevouz with planet point
        end
        
        % Calculation of dv of deep space manoeuvre
        dv=[dv;norm(vds2-vds1)];
        % Initialises variables for next flyby
        time=time+tof; % Time of the next flyby (during istantaneous flyby)
        % vin (velocity vector at the next planet encounter) already initialised
        planet1=planet2; % Next planet
        xp1=xp2; % Position of the next planet
        vp1=vp2; % Velocity of the next planet
    end
end

if options(4) % Plot flag
    axis equal;
    hold off;
end