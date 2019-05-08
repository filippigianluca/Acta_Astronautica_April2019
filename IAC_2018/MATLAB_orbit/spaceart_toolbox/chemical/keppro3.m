function [out,error]=keppro3(in,dt,mu,tol,maxiter)
% Analytic orbit propagator for Keplerian orbits. Uses an algorithm from:
% D. A. Vallado, "Fundamentals of Astrodynamics and Applications, Second
% Edition", Microcosm Press, pp. 101-102.
%
%   [out, error] = keppro3(in, dt, mu [, tol [,maxiter]])
%
% All dimensions to be consistent one another.
%
% INPUT
%   in          Initial state, cartesian coordinates (position, velocity).
%   dt          Propagation time.
%   mu          Planetary constant.
%   tol         Tolerance on the time law. Uses the same number to check if
%               the parameter alpha is zero. If omitted, uses 1e-6.
%   maxiter     Maximum number of itarations for the loop. If omitted, uses
%               1000.
%
% OUTPUT
%   out         Final state, cartesian coordinates (position, velocity).
%   error       Error flag:
%                   0 = no error;
%                   1 = check on (f*gd - fd*g == 1) failed;
%                   2 = number of iterations exceeded before convergence.
%
% FUNCTIONS CALLED
%   (none)
%
% Matteo Ceriotti, 13-02-2007
% Modified by: Matteo Ceriotti, 19/06/2008: Added a check for NaN or Inf.

if dt == 0
    out = in;
    error = 0;
    return
end

if nargin<6
    maxiter=1000;
    if nargin<5
        tol=1e-6;
    end
end

r0=in(1:3);
v0=in(4:6);

nv0=sqrt(v0(1)^2+v0(2)^2+v0(3)^2);
nr0=sqrt(r0(1)^2+r0(2)^2+r0(3)^2);
alpha=-nv0^2/mu+2/nr0;

r0v0=r0(1)*v0(1)+r0(2)*v0(2)+r0(3)*v0(3); % dot(r0, v0)

if alpha>tol % Circle or ellipse
    chi0=sqrt(mu)*dt*alpha;
    if alpha==1
        warning('keppro3:alpha','alpha == 1.');
    end
elseif abs(alpha)<tol % Parabola
    h=[r0(2)*v0(3)-r0(3)*v0(2),r0(3)*v0(1)-r0(1)*v0(3),r0(1)*v0(2)-r0(2)*v0(1)];
    nh=sqrt(h(1)^2+h(2)^2+h(3)^2);
    p=nh^2/mu;
    s=.5*atan(1/(3*sqrt(mu/p^3)*dt));
    w=atan(tan(s)^(1/3));
    chi0=sqrt(p)*2/(tan(2*w));
elseif alpha<-tol
    a=1/alpha;
    chi0=sign(dt)*sqrt(-a)*log(-2*mu*alpha*dt/(r0v0+sign(dt)*sqrt(-mu*a)*(1-nr0*alpha)));
else
    keyboard
end

i=0;
chin=chi0;

% The following lines are to use the tolerance on the time
numer=tol+1;
nr=1;
while abs(numer/nr)>tol && i<=maxiter
% **********
%while abs(chin-chin1)>tol
    i=i+1;
    psi=chin^2*alpha;
    
    if psi>tol
        c2=(1-cos(sqrt(psi)))/psi;
        c3=(sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
    elseif psi<-tol
        c2=(1-cosh(sqrt(-psi)))/psi;
        c3=(sinh(sqrt(-psi))-sqrt(-psi))/sqrt(-psi^3);
    else
        c2=.5;
        c3=1/6;
    end
    %fprintf('%f %f %.12f %.12f\n',chin,psi,c2,c3);
    r=chin^2*c2+r0v0/sqrt(mu)*chin*(1-psi*c3)+nr0*(1-psi*c2);
    if isnan(r) || isinf(r)
        out = [0 0 0 0 0 0];
        error = 3;
        return
    end
    nr=norm(r);
    numer=sqrt(mu)*dt-chin^3*c3-r0v0/sqrt(mu)*chin^2*c2-nr0*chin*(1-psi*c3);
    chin=chin+numer/nr;
end

if i>maxiter
    out = [0 0 0 0 0 0];
    error = 2;
    return
end

f=1-chin^2/nr0*c2;
g=dt-chin^3/sqrt(mu)*c3;
r=f*r0+g*v0;
nr=norm(r);
gd=1-chin^2/nr*c2;
fd=sqrt(mu)/nr/nr0*chin*(psi*c3-1);
v=fd*r0+gd*v0;

out=[r(1),r(2),r(3),v(1),v(2),v(3)];

% Error check
if abs(f*gd-fd*g-1)>.01
    %f*gd-fd*g
    error=1;
else
    error=0;
end

return 