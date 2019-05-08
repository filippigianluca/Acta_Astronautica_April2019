function [Da,De,DI,eDw,DP1,DP2,Dw,DQ1, DQ2]=TD88_Analytic(t,t0,date, X,m,S,solardata, drag, constant)
% function [Da,De,DI,eDw,DP1,DP2,Dw,Draan]=ChangeInOrbitalElements(t,X,m,S,solardata,constant)

% Calculates the change in the orbital elements due to the perturbation
% form the solar activity and geomagnetic sources.
% Input t must be defined in seconds
%       the inital time must be the 1/1/2000 (due to alpha calculation)
% m is the mass of the satellite 
% S is the cross section of the satellite
% Alpha is solar right ascension
% E is the osculating eccentric anomaly

% orbital elements
a=X(1);             % semimajor axis
e=X(2);             % eccentricity
I=X(3);             % inclination
w=X(4);             % argurment of perigee
CapitalOmega=X(5);  % RAAN
f=X(6);             % true anomaly

E=asin((sin(f) * sqrt(1 - (e^2))) / (1 + e * cos(f))); % Valldo page 46


% day count
% d=daycount(t);

current_time = t0 + t * constant.TU / 86400;

% day count in year
% keyboard
d = current_time - date2mjd2000([date.year 1 1 0 0 0]);

while d > 365
    d = d - 365;
end



% Solar right ascension (https://en.wikipedia.org/wiki/Position_of_the_Sun)
% u=juliandate(datetime(t))-juliandate(datetime('01-01-2000','InputFormat','dd-MM-yyyy')); % Modified Julian date | using datetime
% uu=(t - 2000) / constant.one_day;
% epsilon=(23.439 + 0.0000004 * uu) * pi/180;
% LL=(280.46 + 0.9856474 * uu) * pi/180;                     % in radians
% gg=(357.528 + 0.9856003 * uu) * pi/180;                 % in radians
% lambda = LL + (1.915 * pi/180) * sin(gg) + (0.020 * pi/180) * sin(2 * gg);
% Alpha = atan(cos(epsilon) * tan(lambda));

% Alpha = 0;

Alpha= SunRightAscension(current_time);

% local solar time
% ts=localsolartime(X,E,Alpha,t,constant.one_day);
% ts = 1;
% ts = atan2(-cos(I) * sin)



%% 


%% Call scripts 
% For altitudes between 150km and 750km

if  (a-constant.DU * 1e3) <= 150000  && (a-constant.DU * 1e3) > 100000
    [H,c1,p,knj,kn]=TD88_Constants100to150();
elseif (a-constant.DU * 1e3) <= 750000   && (a-constant.DU * 1e3) > 150000 
    [H,c1,p,knj,kn]=TD88_Constants150to750();
else
    [H,c1,p,knj,kn]=TD88_Constants750to1200();
end

% Call Parameters
[fo,fx,ko,g1_4,gB,delta,Tau,cj,zj,Bj]=TD88_Parameters2(X,solardata,drag, constant, d, date,c1,p,H);




% Call derivation Scripts
[Ua,Uja,V5a,V6a,V7a,V5ja,V6ja,V7ja]=Deriv_Da(e,I,w,CapitalOmega,Alpha,Tau,cj,zj,p);

[Ue,Uje,V5e,V6e,V7e,V5je,V6je,V7je]=Deriv_De_moreCorrections(e,I,w,CapitalOmega,Alpha,Tau,cj,zj,p);
[Ui,Uji,V5i,V6i,V7i,V5ji,V6ji,V7ji]=Deriv_DI_corrected(e,I,w,CapitalOmega,Alpha,Tau,cj,zj,p);
[Uw,Ujw,V5w,V6w,V7w,V5jw,V6jw,V7jw]=Deriv_DI_corrected(e,I,w,CapitalOmega,Alpha,Tau,cj,zj,p);

%% Sumation loops 
O=zeros(3,1);

% Da Sum1 
A=zeros(1,3); % empty vector
B=zeros(4,1); % empty vector

for n=1:4
    for j=1:3
    A(j+1)=A(j)+Bj(j)*knj(n,j)*Uja(j); % Summation with j=1:3 for each value of n 
    end
    B(n+1)=g1_4(n)*(kn(n)*Ua+A(j+1)); %summation for n=1:4. this should produce an array which can be sumed.
end

DaSum1=sum(B); 

% Da sum2 
Vna=[0,0,0,0,V5a,V6a,V7a]';
Vnja=[O,O,O,O,V5ja,V6ja,V7ja]';

C=zeros(1,3); % empty vector
D=zeros(4,1); % empty vector

for nn=5:7
    for jj=1:3
    C(jj+1)=C(jj)+Bj(jj)*knj(nn,jj)*Vnja(nn,jj); 
    end
    D(nn+1)=gB(nn)*(kn(nn)*Vna(nn)+C(jj+1)); 
end

DaSum2=sum(D); 



% De Sum1 
A1=zeros(1,3); % empty vector
B1=zeros(4,1); % empty vector

for n=1:4
    for j=1:3
    A1(j+1)=A1(j)+Bj(j)*knj(n,j)*Uje(j); % Summation with j=1:3 for each value of n 
    end
    B1(n+1)=g1_4(n)*(kn(n)*Ue+A1(j+1)); %summation for n=1:4. this should produce an array which can be sumed.
end

DeSum1=sum(B1); 

% De sum2 
Vne=[0,0,0,0,V5e,V6e,V7e]';
Vnje=[O,O,O,O,V5je,V6je,V7je]';

C1=zeros(1,3); % empty vector
D1=zeros(4,1); % empty vector

for nn=5:7
    for jj=1:3
    C1(jj+1)=C1(jj)+Bj(jj)*knj(nn,jj)*Vnje(nn,jj); 
    end
    D1(nn+1)=gB(nn)*(kn(nn)*Vne(nn)+C1(jj+1)); 
end

DeSum2=sum(D1); 




% DI Sum1 
A2=zeros(1,3); % empty vector
B2=zeros(4,1); % empty vector

for n=1:4
    for j=1:3
    A2(j+1)=A2(j)+Bj(j)*knj(n,j)*Uji(j); % Summation with j=1:3 for each value of n 
    end
    B2(n+1)=g1_4(n)*(kn(n)*Ui+A2(j+1)); %summation for n=1:4. this should produce an array which can be sumed.
end

DISum1=sum(B2); 

% DI sum2 
Vni=[0,0,0,0,V5i,V6i,V7i]';
Vnji=[O,O,O,O,V5ji,V6ji,V7ji]';

C2=zeros(1,3); % empty vector
D2=zeros(4,1); % empty vector

for nn=5:7
    for jj=1:3
    C2(jj+1)=C2(jj)+Bj(jj)*knj(nn,jj)*Vnji(nn,jj); 
    end
    D2(nn+1)=gB(nn)*(kn(nn)*Vni(nn)+C2(jj+1)); 
end

DISum2=sum(D2); 




% eDw Sum1 
A3=zeros(1,3); % empty vector
B3=zeros(4,1); % empty vector

for n=1:4
    for j=1:3
    A3(j+1)=A3(j)+Bj(j)*knj(n,j)*Ujw(j); % Summation with j=1:3 for each value of n 
    end
    B3(n+1)=g1_4(n)*(kn(n)*Uw+A3(j+1)); %summation for n=1:4. this should produce an array which can be sumed.
end

eDwSum1=sum(B3); 

% DI sum2
Vnw=[0,0,0,0,V5w,V6w,V7w]';
Vnjw=[O,O,O,O,V5jw,V6jw,V7jw]';

C3=zeros(1,3); % empty vector
D3=zeros(4,1); % empty vector

for nn=5:7
    for jj=1:3
    C3(jj+1)=C3(jj)+Bj(jj)*knj(nn,jj)*Vnjw(nn,jj); 
    end
    D3(nn+1)=gB(nn)*(kn(nn)*Vnw(nn)+C3(jj+1)); 
end

eDwSum2=sum(D3); 

  
%% Period and numerical factor calculation

% P=2*pi/(sqrt(constant.mu/(a^3)));
% TimeStep=86400; % one day in Seconds
% NF = TimeStep/P;
%% Keplarian elements Da, De, DI, eDw

Da=((-1) * 2 * pi * (a^2) * delta * fx * fo * ko * (DaSum1 + DaSum2));
De=((-1) * 2 * pi * a * delta * fx * fo * ko * (DeSum1 + DeSum2));
DI=((-1) * 2 * pi * a * delta * fx * fo * ko * (Tau * tan(I) / (2 * (1 - e^2))) * (DISum1 + DISum2)) ;
eDw=((-1) * 2 * pi * a * delta * fx * fo * ko * (eDwSum1 + eDwSum2)) ;
% Draan=(- 3/2 * sqrt(constant.mu/a^3) * constant.J2 * (constant.Rearth/a)^2 *  cos(I)) * TimeStep;

%% Equinoctical elements 

% % Da
DP1=sin(w+CapitalOmega)*De+cos(w+CapitalOmega)*eDw;
DP2=cos(w+CapitalOmega)*De+sin(w+CapitalOmega)*eDw;


%%% Marilena:
DP1=sin(w)*De+cos(w)*eDw;
DP2=cos(w)*De+sin(w)*eDw;

P1 = e * sin(w);
P2 = e * cos(w);

Dw = cos(w)^2 * (DP1/P2 - P1/P2^2 * DP2);
%%% end Marilena



% DP1=sin(w)*De+cos(w)*eDw;
% DP2=cos(w)*De+sin(w)*eDw;
DQ1=0.5*sin(CapitalOmega)*DI/cos((I/2)^2);
DQ2=0.5*cos(CapitalOmega)*DI/cos((I/2)^2);


end