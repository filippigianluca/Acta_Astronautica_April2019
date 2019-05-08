%--------------------------------------------------------------------------
function [v1,v2,a,p,theta,iter]=lambertI(r1,r2,t,mu,lw,N,branch)
%
%This routine implements a new algorithm that solves Lambert's problem. The
%algorithm has two major characteristics that makes it favorable to other
%existing ones.
%
%   1) It describes the generic orbit solution of the boundary condition
%   problem through the variable X=log(1+cos(alpha/2)). By doing so the
%   graphs of the time of flight become defined in the entire real axis and
%   resembles a straight line. Convergence is granted within few iterations
%   for all the possible geometries (except, of course, when the transfer
%   angle is zero). When multiple revolutions are considered the variable is
%   X=tan(cos(alpha/2)*pi/2).
%
%   2) Once the orbit has been determined in the plane, this routine
%   evaluates the velocity vectors at the two points in a way that is not
%   singular for the transfer angle approaching to pi (Lagrange coefficient
%   based methods are numerically not well suited for this purpose).
%
%   As a result Lambert's problem is solved (with multiple revolutions
%   being accounted for) with the same computational effort for all
%   possible geometries. The case of near 180 transfers is also solved
%   efficiently.
%
%   We note here that even when the transfer angle is exactly equal to pi
%   the algorithm does solve the problem in the plane (it finds X), but it
%   is not able to evaluate the plane in which the orbit lies. A solution
%   to this would be to provide the direction of the plane containing the
%   transfer orbit from outside. This has not been implemented in this
%   routine since such a direction would depend on which application the
%   transfer is going to be used in.
%
%Usage: [v1,v2,a,p,theta,iter]=lambertI(r1,r2,t,mu,lw,N,branch)
%
%Inputs:
%           r1=Position vector at departure (column)
%           r2=Position vector at arrival (column, same units as r1)
%           t=Transfer time (scalar)
%           mu=gravitational parameter (scalar, units have to be
%           consistent with r1,t units)
%           lw=1 if long way is chosen
%           branch='l' if the left branch is selected in a problem where N
%           is not 0 (multirevolution)
%           N=number of revolutions
%
%Outputs:
%           v1=Velocity at departure        (consistent units)
%           v2=Velocity at arrival
%           a=semi major axis of the solution
%           p=semi latus rectum of the solution
%           theta=transfer angle in rad
%           iter=number of iteration made by the newton solver (usually 6)
%
%please report bugs to dario.izzo@esa.int





%Preliminary control on the function call
if nargin==5
    N=0;
end
if t<=0
    warning('Negative time as input')
    v1=NaN;
    v2=NaN;
    return
end


tol=1e-11;  %Increasing the tolerance does not bring any advantage as the
%precision is usually greater anyway (due to the rectification of the tof
%graph) except near particular cases such as parabolas in which cases a
%lower precision allow for usual convergence.


%Non dimensional units
R=sqrt(r1'*r1);
V=sqrt(mu/R);
T=R/V;

%working with non-dimensional radii and time-of-flight
r1=r1/R;
r2=r2/R;
t=t/T;

%Evaluation of the relevant geometry parameters in non dimensional units
r2mod=sqrt(r2'*r2);
theta=real(acos((r1'*r2)/r2mod)); %the real command is useful when theta is very
%close to pi and the acos function could return complex numbers
if lw
    theta=2*pi-theta;
end
c=sqrt(1+r2mod^2-2*r2mod*cos(theta)); %non dimensional chord
s=(1+r2mod+c)/2;                      %non dimensional semi-perimeter
am=s/2;                               %minimum energy ellipse semi major axis
lambda=sqrt(r2mod)*cos(theta/2)/s;    %lambda parameter defined in BATTIN's book



%We start finding the log(x+1) value of the solution conic:
%%NO MULTI REV --> (1 SOL)
if N==0
    inn1=-.5233;    %first guess point
    inn2=.5233;     %second guess point
    x1=log(1+inn1);
    x2=log(1+inn2);
    y1=log(x2tof(inn1,s,c,lw,N))-log(t);
    y2=log(x2tof(inn2,s,c,lw,N))-log(t);

    %Newton iterations
    err=1;
    i=0;
    while ((err>tol) && (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=log(x2tof(exp(xnew)-1,s,c,lw,N))-log(t);
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);
    end
    iter=i;
    x=exp(xnew)-1;


    %%MULTI REV --> (2 SOL) SEPARATING RIGHT AND LEFT BRANCH
else
    if branch=='l'
        inn1=-.5234;
        inn2=-.2234;
    else
        inn1=.7234;
        inn2=.5234;
    end
    x1=tan(inn1*pi/2);
    x2=tan(inn2*pi/2);
    y1=x2tof(inn1,s,c,lw,N)-t;

    y2=x2tof(inn2,s,c,lw,N)-t;
    err=1;
    i=0;

    %Newton Iteration
    while ((err>tol) && (i<60) && (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=x2tof(atan(xnew)*2/pi,s,c,lw,N)-t;
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);
    end
    x=atan(xnew)*2/pi;
    iter=i;
end

%The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
%now need the conic. As for transfer angles near to pi the lagrange
%coefficient technique goes singular (dg approaches a zero/zero that is
%numerically bad) we here use a different technique for those cases. When
%the transfer angle is exactly equal to pi, then the ih unit vector is not
%determined. The remaining equations, though, are still valid.


a=am/(1-x^2);                       %solution semimajor axis
%calcolo psi
if x<1 %ellisse
    beta=2*asin(sqrt((s-c)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acos(x);
    psi=(alfa-beta)/2;
    eta2=2*a*sin(psi)^2/s;
    eta=sqrt(eta2);
else %iperbole
    beta=2*asinh(sqrt((c-s)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acosh(x);
    psi=(alfa-beta)/2;
    eta2=-2*a*sinh(psi)^2/s;
    eta=sqrt(eta2);
end
p=r2mod/am/eta2*sin(theta/2)^2;     %parameter of the solution
sigma1=1/eta/sqrt(am)*(2*lambda*am-(lambda+x*eta));
ih=vers(vett(r1,r2)');
if lw
    ih=-ih;
end

vr1 = sigma1;
vt1 = sqrt(p);
v1  = vr1 * r1   +   vt1 * vett(ih,r1)';

vt2=vt1/r2mod;
vr2=-vr1+(vt1-vt2)/tan(theta/2);
v2=vr2*r2/r2mod+vt2*vett(ih,r2/r2mod)';
v1=v1*V;
v2=v2*V;
a=a*R;
p=p*R;


%--------------------------------------------------------------------------
function t=x2tof(x,s,c,lw,N)
%Subfunction that evaluates the time of flight as a function of x
am=s/2;
a=am/(1-x^2);
if x<1 %ELLISSE
    beta=2*asin(sqrt((s-c)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acos(x);
else   %IPERBOLE
    alfa=2*acosh(x);
    beta=2*asinh(sqrt((s-c)/(-2*a)));
    if lw
        beta=-beta;
    end
end
t=tofabn(a,alfa,beta,N);



%--------------------------------------------------------------------------
function t=tofabn(sigma,alfa,beta,N)
%subfunction that evaluates the time of flight via Lagrange expression
if sigma>0
    t=sigma*sqrt(sigma)*((alfa-sin(alfa))-(beta-sin(beta))+N*2*pi);
else
    t=-sigma*sqrt(-sigma)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
end


%--------------------------------------------------------------------------
function v=vers(V)
%subfunction that evaluates unit vectors
v=V/sqrt(V'*V);

%--------------------------------------------------------------------------
function [DV,rp,iter]=PowSwingByInv(Vin,Vout,alpha)
%
%Usage: [DV,rp,iter] = PowSwingByInv(Vin,Vout,alpha)
%
%Outputs:
%           DV:    Velcity Increment of the Powered SwingBy (non dimensional)
%           rp:    Pericenter radius found.
%           iter:  Number of iteration to converge (-1 if convergence is failed)
%
%Inputs:
%           Vin:   Incoming hyperbolic velocity modulus  (non dimensional)
%           Vout:  Outgoing hyperbolic velocity modulus  (non dimensional)
%           alpha: Angle between Vin and Vout (in rad.)
%
%Comments:  The non dimensional units are R for the length and sqrt(mu/R)
%for the velocity --> gravitational constant is one. R may be choosen
%freely as any relevant length. Magic of the non dimensional forms: if we
%forget about dimension and we call the routine, then DV is returned in the
%same units as the input parameters, and rp has to be multiplied by the
%planet gravitational constant (unit consistent with the velocity input)
%to be transformed in the length.

aIN=1/Vin^2;    %semimajor axis of the incoming hyperbola
aOUT=1/Vout^2;  %semimajor axis of the outcoming hyperbola

%We find the perigee radius with an iteration method based on the gradient
%of the function. Attention has to be payed to the initial point as the
%function is not defined for rp<0. The option here implemented considers
%halfing the perigee radius whenever the gradient pushes the next iteration
%in a non defined zone.
i=0;
maxiter=30;     %maximum number of iteration allowed for the gradient method
rp=1;           %Initial point
err=1;
while((err>1e-8)&&(i<maxiter))
    i=i+1;
    f=asin(aIN/(aIN+rp))+asin(aOUT/(aOUT+rp))-alpha;
    df=-aIN/sqrt(rp^2+2*aIN*rp)/(aIN+rp)-aOUT/sqrt(rp^2+2*aOUT*rp)/(aOUT+rp);
    rpNew=rp-f/df;
    if (rpNew>0)
        err=abs(rpNew-rp);
        rp=rpNew;
    else
        rp=rp/2;
    end
end

%Evaluation of the DV
DV=abs(sqrt(Vout^2+2/rp)-sqrt(Vin^2+2/rp));

%If the maximum number of iteration is achieved the returned number of
%iteration is -1.
iter=i;
if iter==maxiter
    iter=-1;
end


%--------------------------------------------------------------------------
function [r,v]=CUSTOMeph(jd,epoch,keplerian,flag)
%
%Returns the position and the velocity of an object having keplerian
%parameters epoch,a,e,i,W,w,M
%
%Usage:     [r,v]=CUSTOMeph(jd,name,list,data,flag)
%
%Inputs:    jd: julian date
%           epoch: mjd when the object was observed (referred to M)
%           keplerian: vector containing the keplerian orbital parameters
%
%Output:    r = object position with respect to the Sun (km if flag=1, AU otherwise)
%           v = object velocity ( km/s if flag=1, AU/days otherways )
%
%Revisions :    Function added 04/07

muSUN=1.32712428e+11;    %Gravitational constant of Sun
AU  = 149597870.66;      %Astronomical Unit in Km



     a=keplerian(1)*AU; %in km
     e=keplerian(2);
     i=keplerian(3); 
     W=keplerian(4);
     w=keplerian(5);
     M=keplerian(6);
     jdepoch=mjd2jed(epoch);
     DT=(jd-jdepoch)*60*60*24;
     n=sqrt(muSUN/a^3);
     M=M/180*pi;
     M=M+n*DT;
     M=mod(M,2*pi);
     E=M2E(M,e);
     [r,v]=par2IC([a,e,i/180*pi,W/180*pi,w/180*pi,E],muSUN);
     if flag~=1
         r=r/AU;
         v=v*86400/AU;
     end
 

%--------------------------------------------------------------------------
function jd = mjd20002jed(mjd2000)
%This function converts mean julian date 2000 to julian date
jd=mjd2000+2451544.5;

%--------------------------------------------------------------------------
function jd = mjd2jed(mjd)
% This function converts mean Julian date into Julian date
jd = mjd +2400000.5;

%--------------------------------------------------------------------------
function [r0,v0]=par2IC(E,mu)
%
%Usage: [r0,v0] = IC2par(E,mu)
%
%Outputs:
%           r0:    column vector for the position
%           v0:    column vector for the velocity
%
%Inputs:
%           E:     Column Vectors containing the six keplerian parameters,
%                  (a (negative fr hyperbolas),e,i,OM,om,Eccentric Anomaly 
%                    or Gudermannian if e>1)
%           mu:    gravitational constant
%
%Comments:  The parameters returned are, of course, referred to the same
%ref. frame in which r0,v0 are given. a can be given either in kms or AUs,
%but has to be consistent with mu.All the angles must be given in radians. 

a=E(1);
e=E(2);
i=E(3);
omg=E(4);
omp=E(5);
EA=E(6);


% Grandezze definite nel piano dell'orbita
if e<1
    b=a*sqrt(1-e^2);
    n=sqrt(mu/a^3);

    xper=a*(cos(EA)-e);
    yper=b*sin(EA);

    xdotper=-(a*n*sin(EA))/(1-e*cos(EA));
    ydotper=(b*n*cos(EA))/(1-e*cos(EA));
else
    b=-a*sqrt(e^2-1);
    n=sqrt(-mu/a^3);   
    dNdzeta=e*(1+tan(EA)^2)-(1/2+1/2*tan(1/2*EA+1/4*pi)^2)/tan(1/2*EA+1/4*pi);
    
    xper = a/cos(EA)-a*e;
    yper = b*tan(EA);
    
    xdotper = a*tan(EA)/cos(EA)*n/dNdzeta;
    ydotper = b/cos(EA)^2*n/dNdzeta;
end

% Matrice di trasformazione da perifocale a ECI

R(1,1)=cos(omg)*cos(omp)-sin(omg)*sin(omp)*cos(i);
R(1,2)=-cos(omg)*sin(omp)-sin(omg)*cos(omp)*cos(i);
R(1,3)=sin(omg)*sin(i);
R(2,1)=sin(omg)*cos(omp)+cos(omg)*sin(omp)*cos(i);
R(2,2)=-sin(omg)*sin(omp)+cos(omg)*cos(omp)*cos(i);
R(2,3)=-cos(omg)*sin(i);
R(3,1)=sin(omp)*sin(i);
R(3,2)=cos(omp)*sin(i);
R(3,3)=cos(i);

% Posizione nel sistema inerziale 

r0=R*[xper;yper;0];
v0=R*[xdotper;ydotper;0];