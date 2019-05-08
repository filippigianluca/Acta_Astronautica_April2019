function [History,F_n,Inclination,dF_n,di_dt,RAAN] = BinormalControl_analytical(DepElements,ArrElements,r,G,range,p,nSteps_i_h)

%% Define useful parameters
    a0= p(1); a1= p(2); a2= p(3); a3 = p(4); a4 = p(5); a5 = p(6); a6=p(7);
    G0 = p(8); G1 = p(9); G2 = p(10); G3 = p(11);

u_v = [DepElements(4) ArrElements(4)];   

    
%% Number of line of node passages (For inclination change)
i_0 = DepElements(5);
i_f = ArrElements(5);

InitialIncStep = i_0;
DeltaInclinationToFinalValue = i_f - i_0;

%% Number of u = 90-270 passages (For inclination change)
h_0 = DepElements(6);
h_f = ArrElements(6);

InitialHStep = h_0 ;
DeltaHToFinalValue = h_f - h_0;

%% Find where i and h do not change

[U,i_or_h] = TrajectoryInspection (u_v(1),u_v(end),range);

% Number of steps
TotalUHChange = (U(:,2)-U(:,1))'*i_or_h;

TotalUInclinationChange = (U(:,2)-U(:,1))'*(i_or_h == 0);

%% Compute F_n

j = 0; % For filling the RESULT matrix 

F_n = @(u) 0 ;
Inclination = @(u) 0 ;
RAAN = @(u) 0;
di_dt = @(u) 0 ;
dF_n = @(u) 0;
History = [];

if U(1,1) > u_v(1)+1e-3
    
    j = j+1 ;

%     RESULT{j,1} = i_0;
%     RESULT{j,2} = h_0;
%     RESULT{j,3} = 0;
%     RESULT{j,4} = [u_v(1) U(1,1)];
%     RESULT{j,5} = 0;
%     
    History = [u_v(1)];

    Inclination = @(u) Inclination(u) + i_0.*( u>=u_v(1) & u<U(1,1) ) ;
    RAAN = @(u) RAAN(u)+h_0.*( u>=u_v(1) & u<U(1,1) ) ;
end

ll = size(U,1);

for ii = 1:ll
    
    if i_or_h(ii) == 0

        dU = U(ii,2) - U(ii,1);
        
        NumberInclinationSteps = dU ./ TotalUInclinationChange;  
        
        DeltaInclinationInThisStep = DeltaInclinationToFinalValue .* NumberInclinationSteps ;
        
        FinalIncStep = InitialIncStep + DeltaInclinationInThisStep;
        
        UIntervalStep = linspace(U(ii,1),U(ii,2),nSteps_i_h)';
        du = UIntervalStep(2) - UIntervalStep(1) ;
        
        % Define boundary conditions
        u0 = UIntervalStep(1);
        i0 = InitialIncStep ;
        i1 = 0;
        i_aux = [dU^3 dU^2 ; 3*dU^2 2*dU]\[FinalIncStep-InitialIncStep ; 0 ] ;
        i2 = i_aux(2) ; 
        i3 = i_aux(1) ;
        InclinationShapeStep = @(u) i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3 ; 
        Inclination = @(u) Inclination(u) + InclinationShapeStep(u).*( u>= U(ii,1) & u<U(ii,2)) ;
        di_dt = @(u) di_dt(u) + ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2 ).*( u>= U(ii,1) & u <= U(ii,2)) ;
        
        
        %% Try with optimal control theory
%         u_2pi = wrapToPiHalf(UIntervalStep) ;
%         u_2pi = UIntervalStep;
%         r3G2 = r(UIntervalStep).^3./G(UIntervalStep).^2;
%         pp = polyfit(u_2pi,r3G2,1);
%         A = pp(1) ; B = pp(2);
%         yy = A.*(u_2pi)+B;
%         
%         max_error = max(abs((r3G2-yy)./r3G2 *100));
%         if max_error > 1
%             keyboard
%         end
%      
% %         alpha0 = wrapToPiHalf(U(ii,1)) ;
%         alpha0 = U(ii,1) ;
% %         alphaf = wrapToPiHalf(U(ii,2)) ;
%         alphaf = U(ii,2) ;
%         
%         a = @(u) 0.5*( u.^2 .* sin(u) .* cos(u) + u .^ 3 + u .* cos(u) .^2) - 0.25.*(cos(u) .* sin(u) + u) - (u.^3)./3 ;
%         c = @(u) 0.5 .* u .* ( sin(u) .* cos(u) + u ) - 0.25 .* ( sin(u) .^ 2 + u .^ 2 ) ;
%         b = @(u) 0.5 .* ( cos(u) .* sin(u) + u );
%         A0 = a(alpha0);
%         C0 = c(alpha0);
%         B0 = b(alpha0);
%         Af = a(alphaf);
%         Cf = c(alphaf);
%         Bf = b(alphaf);
%         lambda_x = 2 * (FinalIncStep-InitialIncStep) / ( A^2 * (A0 - Af) + B^2 * (B0 - Bf) + 2*A*B * (C0 - Cf) ) ;
%         K = InitialIncStep + lambda_x / 2 * (A0*A^2 + B0*B^2 + 2*C0*A*B);
%         K2 = FinalIncStep + lambda_x / 2 * (Af*A^2 + Bf*B^2 + 2*Cf*A*B);
%         
%         if K~=K2
%             keyboard
%         end
%         
%         Inclination_try = @(u) -lambda_x / 2 * (A .^2 .* a(u) + 2.*A.*B .* c(u) + b(u) .* B.^2) + K;
%         F_n_try = @(u) - lambda_x / 2 .* ( A.*u + B ) .* cos(u) ;
%         
%         [t,x] = RK4_4 (UIntervalStep,u_2pi,InitialIncStep,F_n_try,r,G);
        
        
%         i2 = @(u) 
        
        
        % Find Control of current step
%         F_n_Step = @(u) ( i1 .* G(u) .^ 2 ./ r(u) .^ 3 ./ (cos(u) + cot(i0 + i1 .* (u - u0)) .* sin(u) .* i1) );
        F_n_Step = @(u)     (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) ./ (cos(u) + cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .* sin(u) .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2));
        F_n = @(u) F_n(u) + F_n_Step(u) .* (u >= U(ii,1) &  u <= U(ii,2)) ;
%         dF_n = @(u) dF_n(u) + (u == UIntervalStep(1)).*(F_n(U(ii,1))/du)+ (u == UIntervalStep(end)).*(-F_n(UIntervalStep(end-1))/du) + (u > UIntervalStep(1) &  u < U(ii,2)).*(0.2e1 .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) ./ (cos(u) + cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .* sin(u) .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .* (G1 + G2 .* cos(u - u_v(1) + G3)) + 0.3e1 .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) ./ (cos(u) + cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .* sin(u) .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) ./ (cos(u) + cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .* sin(u) .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .^ 2 .* (-sin(u) + ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) .^ 2 .* (-0.1e1 - cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .^ 2) .* sin(u) + cot(i0 + i1 .* (u - u0) + i2 .* (u - u0).^2 + i3 .* (u - u0) .^ 3) .* cos(u) .* ( i1 + 2.* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)));
        dF_n = @(u) dF_n(u) + (u >= U(ii,1) &  u <= U(ii,2)).*( (2 .* i2 + 6 .* i3 .* (u - u0)) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (cos(u) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* sin(u) .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) + 0.2e1 .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (cos(u) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* sin(u) .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .* (G1 + G2 .* cos((u - u_v(1) + G3))) + 0.3e1 .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 ./ (cos(u) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* sin(u) .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .* (a1 + (2 .* a2 .* u) + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (cos(u) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* sin(u) .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2)) .^ 2 .* (-sin(u) + ((i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) .^ 2) .* (-0.1e1 - cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .^ 2) .* sin(u) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* cos(u) .* (i1 + 2 .* i2 .* (u - u0) + 3 .* i3 .* (u - u0) .^ 2) + cot((i0 + i1 .* (u - u0) + i2 .* (u - u0) .^ 2 + i3 .* (u - u0) .^ 3)) .* sin(u) .* (2 .* i2 + 6 .* i3 .* (u - u0))) );
        
        
%         u_dot = G(UIntervalStep)./r(UIntervalStep).^2 - r(UIntervalStep).*cot(Inclination(UIntervalStep)).*sin(UIntervalStep).*F_n_Step(UIntervalStep)./G(UIntervalStep);
%         dvec = abs(F_n_Step(UIntervalStep)./u_dot);
%         dvec = dvec(1:end-1) ;
%         dv1 = (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
%             + 4*sum(dvec(2:2:end-1)) + dvec(end)) .*  149597870.66 ./ 86164.1004
        
%         u_dot = G(UIntervalStep)./r(UIntervalStep).^2 - r(UIntervalStep).*cot(Inclination_try(u_2pi)).*sin(UIntervalStep).*F_n_try(UIntervalStep)./G(UIntervalStep);
%         dvec = abs(F_n_try(UIntervalStep)./u_dot);
%         dvec = dvec(1:end-1) ;
%         dv2 = (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
%             + 4*sum(dvec(2:2:end-1)) + dvec(end)) .*  149597870.66 ./ 86164.1004
%         
%         keyboard
       
        sol = ode45(@(u,x) DiffInclinationHEquationsWithControl(u,x,r,G,F_n_Step,1,InclinationShapeStep),UIntervalStep,InitialHStep);
        [hStep] = deval(sol,UIntervalStep)';
        FinalHStep = hStep(end) ;
        
        pph = polyfit(UIntervalStep-u0,hStep,4);
%         INC = @(u) ppp(1).*(u-u0).^6+ppp(2).*(u-u0).^5+ppp(3).*(u-u0).^4+ppp(4).*(u-u0).^3+ppp(5).*(u-u0).^2+ppp(6).*(u-u0)+ppp(7);
        RAAN = @(u) RAAN(u)+ ( u>= U(ii,1) & u<U(ii,2)).*(pph(1).*(u-u0).^4+pph(2).*(u-u0).^3+pph(3).*(u-u0).^2+pph(4).*(u-u0)+pph(5));
        
        j = j + 1 ; 
        History = [History; UIntervalStep ];

        j = j+1 ;
        
        if U(ii,2) < u_v(end)-1e-6

            if ii == ll
                History = [History ;  u_v(end) ];
                Inclination = @(u) Inclination(u) + FinalIncStep.*( u>=U(ii,2) & u<= u_v(end) ) ;
                RAAN = @(u) RAAN(u) + FinalHStep.*( u>=U(ii,2) & u<= u_v(end) ) ;
            else

                Inclination = @(u) Inclination(u) + FinalIncStep.*( u>=U(ii,2) & u<U(ii+1,1) ) ;
                RAAN = @(u) RAAN(u) + FinalHStep.*( u>=U(ii,2) & u<U(ii+1,1) ) ;
                
            end

        end
        
        InitialIncStep = FinalIncStep ;
        DeltaInclinationToFinalValue = i_f - InitialIncStep ;
         
        InitialHStep = FinalHStep ;
        DeltaHToFinalValue = h_f - InitialHStep;
        
        TotalUInclinationChange = TotalUInclinationChange - dU ;
        
    elseif i_or_h(ii) == 1
        
 
        dU = U(ii,2) - U(ii,1);
        
        NumberHSteps = dU ./ TotalUHChange ;  
        
        DeltaHInThisStep = DeltaHToFinalValue .* NumberHSteps ;
        
        FinalHStep = InitialHStep + DeltaHInThisStep;
        
        UIntervalStep = linspace(U(ii,1),U(ii,2),nSteps_i_h)';
        du = UIntervalStep(2) - UIntervalStep(1);
                    
        % Define boundary conditions
        
        u0 = UIntervalStep(1);
        h0 = InitialHStep ;
        h1 = 0;
        h_aux = [dU^3 dU^2 ; 3*dU^2 2*dU]\[FinalHStep-InitialHStep ; 0] ;
        h2 = h_aux(2) ;
        h3 = h_aux(1) ;
        HShapeStep = @(u) h0 + h1 .* (u - u0) + h2 .* (u - u0) .^ 2 + h3 .* (u - u0) .^ 3 ;
        RAAN = @(u) RAAN(u) + HShapeStep(u).*(u >= U(ii,1) &  u < U(ii,2));
        
        F_n_Step = @(u) (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(InitialIncStep) + cot(InitialIncStep) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2));       
        iStep_old = InitialIncStep .* ones(length(UIntervalStep),1);
        
        %% PUOI PENSARE A UN ITERATIVE METHOD
        
        sol = ode45(@(u,x) DiffInclinationHEquationsWithControl(u,x,r,G,F_n_Step,0),UIntervalStep,InitialIncStep);
        [iStep_new] = deval(sol,UIntervalStep)';
        
        jj = 0 ;

        while max(abs(iStep_new-iStep_old))>1e-12 && jj < 100
        jj = jj +1 ;
        iStep_old = iStep_new;
        ppp = polyfit(UIntervalStep-u0,iStep_old,4);
%         INC = @(u) ppp(1).*(u-u0).^6+ppp(2).*(u-u0).^5+ppp(3).*(u-u0).^4+ppp(4).*(u-u0).^3+ppp(5).*(u-u0).^2+ppp(6).*(u-u0)+ppp(7);
        INC = @(u) ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5);
%   Interpolo con degree 4  
        F_n_Step = @(u) (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(INC(u)) + cot(INC(u)) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2));       
        
        sol = ode45(@(u,x) DiffInclinationHEquationsWithControl(u,x,r,G,F_n_Step,0),UIntervalStep,InitialIncStep);
        [iStep_new] = deval(sol,UIntervalStep)';
        end
       
        if jj == 100
            error('Non interpolato bene')
        end
        FinalIncStep = iStep_new(end) ;
        F_n = @(u) F_n(u) + F_n_Step(u) .* (u >= U(ii,1) &  u <= U(ii,2)) ;
%         dF_n = @(u) dF_n(u) + (u == UIntervalStep(1)).*(F_n(U(ii,1))/du) + (u == UIntervalStep(end)).*(-F_n(UIntervalStep(end-1))/du) + (u > UIntervalStep(1) &  u < U(ii,2)).*(0.2e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) + cot(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (G1 + G2 .* cos(u - u_v(1) + G3)) + 0.3e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 ./ (sin(u) ./ sin(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) + cot(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_v(1)) + G2 .* sin(u - u_v(1) + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) + cot(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .^ 2 .* (cos(u) ./ sin(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) - sin(u) ./ sin(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .^ 2 .* (0.6e1 .* ppp(1) .* (u - u0) .^ 5 + 0.5e1 .* ppp(2) .* (u - u0) .^ 4 + 0.4e1 .* ppp(3) .* (u - u0) .^ 3 + 0.3e1 .* ppp(4) .* (u - u0) .^ 2 + 0.2e1 .* ppp(5) .* (u - u0) + ppp(6)) .* cos(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) + (0.6e1 .* ppp(1) .* (u - u0) .^ 5 + 0.5e1 .* ppp(2) .* (u - u0) .^ 4 + 0.4e1 .* ppp(3) .* (u - u0) .^ 3 + 0.3e1 .* ppp(4) .* (u - u0) .^ 2 + 0.2e1 .* ppp(5) .* (u - u0) + ppp(6)) .* (-0.1e1 - cot(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .^ 2) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) + cot(ppp(1) .* (u - u0) .^ 6 + ppp(2) .* (u - u0) .^ 5 + ppp(3) .* (u - u0) .^ 4 + ppp(4) .* (u - u0) .^ 3 + ppp(5) .* (u - u0) .^ 2 + ppp(6) .* (u - u0) + ppp(7)) .* cos(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)));
   
        % Interpola con degree 6
%         dF_n = @(u) dF_n(u) + (u > U(ii,1) &  u <= U(ii,2)).*( (2 .* h2 + 6 .* h3 .* (u - u0)) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) + 0.2e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (G1 + G2 .* cos((u - u_v(1) + G3))) + 0.3e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (a1 + (2 .* a2 .* u) + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .^ 2 .* (cos(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) - sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .^ 2 .* ((0.6e1 .* ppp(1) .* (u - u0) .^ 5 + 0.5e1 .* ppp(2) .* (u - u0) .^ 4 + 0.4e1 .* ppp(3) .* (u - u0) .^ 3 + 0.3e1 .* ppp(4) .* (u - u0) .^ 2 + 0.2e1 .* ppp(5) .* (u - u0) + ppp(6))) .* cos((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + ((0.6e1 .* ppp(1) .* (u - u0) .^ 5 + 0.5e1 .* ppp(2) .* (u - u0) .^ 4 + 0.4e1 .* ppp(3) .* (u - u0) .^ 3 + 0.3e1 .* ppp(4) .* (u - u0) .^ 2 + 0.2e1 .* ppp(5) .* (u - u0) + ppp(6))) .* (-0.1e1 - cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .^ 2) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* cos(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (2 .* h2 + 6 .* h3 .* (u - u0))) ) ;

        %   Interpola INC con degree 4      
        dF_n = @(u) dF_n(u) + (u >= U(ii,1) &  u <= U(ii,2)).*( (2 .* h2 + 6 .* h3 .* (u - u0)) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) + 0.2e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (G1 + G2 .* cos((u - u_v(1) + G3))) + 0.3e1 .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .* (a1 + (2 .* a2 .* u) + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + (G1 .* (u - u_v(1))) + G2 .* sin((u - u_v(1) + G3))) .^ 2 .* (a0 + (a1 .* u) + (a2 .* u .^ 2) + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2)) .^ 2 .* (cos(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) - sin(u) ./ sin((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .^ 2 .* ((4 .* ppp(1) .* (u - u0) .^ 3 + 3 .* ppp(2) .* (u - u0) .^ 2 + 2 .* ppp(3) .* (u - u0)  + ppp(4))) .* cos((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) + ((4 .* ppp(1) .* (u - u0) .^ 3 + 3 .* ppp(2) .* (u - u0) .^ 2 + 2 .* ppp(3) .* (u - u0)  + ppp(4))) .* (-0.1e1 - cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .^ 2) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* cos(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) + cot((ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5))) .* sin(u) .* (2 .* h2 + 6 .* h3 .* (u - u0))) ) ;
         
        Inclination = @(u) Inclination(u) + INC(u).*( u>= U(ii,1) & u < U(ii,2)) ;
        
        % Degree 6
%         di_dt = @(u) di_dt(u) + ( u>= U(ii,1) & u <= U(ii,2)).*(6*ppp(1).*(u-u0).^5+5*ppp(2).*(u-u0).^4+4*ppp(3).*(u-u0).^3+3*ppp(4).*(u-u0).^2+2*ppp(5).*(u-u0)+ppp(6)) ;
        
        % Degree 4
        di_dt = @(u) di_dt(u) + ( u>= U(ii,1) & u <= U(ii,2)).*(4.*ppp(1).*(u-u0).^3+3.*ppp(2).*(u-u0).^2+2.*ppp(3).*(u-u0)+ppp(4)) ;
        
%       [F_n_dis,INC] = F_n_Step_Calculator(h1,InitialIncStep,r,G,[linspace(UIntervalStep(1),UIntervalStep(end),1000000)]);

%       r2d = 180/pi ;
%       plot(linspace(UIntervalStep(1),UIntervalStep(end),1000000)*r2d,INC*r2d,'r',UIntervalStep*r2d,iStep*r2d,'b')
%
%
        j = j+1;
        History = [History; UIntervalStep ];
        j = j+1 ;
        if U(ii,2) < u_v(end)-1e-6

            if ii == ll
                History = [History ;  u_v(end) ];
                Inclination = @(u) Inclination(u) + FinalIncStep.*( u>=U(ii,2) & u<= u_v(end) ) ;
                RAAN = @(u) RAAN(u) + FinalHStep.*( u>=U(ii,2) & u<= u_v(end) ) ;
            else

                Inclination = @(u) Inclination(u) + FinalIncStep.*( u>=U(ii,2) & u<U(ii+1,1) ) ;
                RAAN = @(u) RAAN(u) + FinalHStep.*( u>=U(ii,2) & u<U(ii+1,1) ) ;
            end

        end
        
        InitialIncStep = FinalIncStep ;
        DeltaInclinationToFinalValue = i_f - InitialIncStep ;
        
        
        InitialHStep = FinalHStep ;
        DeltaHToFinalValue = h_f - InitialHStep ;
        
        TotalUHChange = TotalUHChange - dU ;
       
    else
        error('Unknown error')        
    end
    
    
end
    


end
    
