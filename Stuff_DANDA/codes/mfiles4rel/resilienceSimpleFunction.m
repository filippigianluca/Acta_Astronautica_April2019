function RES = resilienceSimpleFunction(TM,r2,r1,lam, mu)
%TM is the mission time
%r2 is the immediate performance function of fully operational system
%r1 is the immediate performance function of partially failed system
%lam is the rate of transition to the partially failed state
%mu is the rate of transition to the fully operational state

%this has been replaced by ``satellite_survival'' function
%S is the survival function of the satelite (probability that it has not failed completely yet)
% e.g. S = @(x) 1-wblcdf(x,5000,0.8)

function prob = p2(t)
    s = mu+lam;
    prob = mu/s + lam*exp(-t*s)/s;
end

function R = ER(t)
    pff = p2(t);
%     R = ( r2(t).*pff +r1(t).*(1.-pff) ).*satellite_survival(t);
    R = ( r2.*pff +r1.*(1.-pff) ).*satellite_survival(t);
end

% t = 0:TM/150:TM;
% y = ER(t);
% plot(t,y);

RES = quad(@ER,0,TM, 0.3);%/TM;

% RES = integral(@ER,0,TM);%/TM;
%     function rr = r2fun(t)
%         rr = r2(t);
%     end
end