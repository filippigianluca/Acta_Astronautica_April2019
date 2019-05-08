function [DV,TIME] = BinormalControl(DepElements,ArrElements,range,p,mu,du_i_h,du_others)

%% Define useful parameters
a0= p(1); a1= p(2); a2= p(3); a3 = p(4); a4 = p(5); a5 = p(6); a6=p(7);
G0 = p(8); G1 = p(9); G2 = p(10); G3 = p(11);

u_0 = DepElements(4); u_f = ArrElements(4) ;

r = @(u) 1./( a0 + a1.*u + a2.*u.^2 + ( a3 + a4.*u ).*cos(u) + ( a5 + a6.*u ).*sin(u) );

G = @(u) G0 + G1.*(u-u_0) + G2.*sin( u-u_0 + G3 );

% 
% U = linspace(u_0,u_f,100000);
% figure
% plot(U,r(U))
% grid on
% 
% figure
% plot(U,G(U))
% grid on
%     
% 
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
range = abs(range) ;
[U,i_or_h] = TrajectoryInspection (u_0,u_f,range);

ll = size(U,1);

% Number of steps
TotalUHChange = (U(:,2)-U(:,1))'*i_or_h;

TotalUInclinationChange = (U(:,2)-U(:,1))'*(i_or_h == 0);

%% Compute F_n

DV = 0;
TIME = 0;

if U(1,1) > u_0+1e-3
    
    n_steps = ceil((U(1,1)-du_others - u_0)/du_others +1) ;
    if n_steps < 5
        n_steps = 5 ;
    end
    UIntervalStep = linspace(u_0,U(1,1),n_steps+1)';
    UIntervalStep = UIntervalStep(1:end-1);
    du = UIntervalStep(2)-UIntervalStep(1);
    
    sin_u = sin(UIntervalStep);
    cos_u = cos(UIntervalStep);
            
    r_step = r(UIntervalStep) ;
    dr_step = - r_step .^2 .* ( a1 + 2 .* a2 .* UIntervalStep + a4 .* cos_u - ( a4 .* UIntervalStep + a3 ) .* sin_u + a6 .* sin_u + ( a6 .* UIntervalStep + a5 ) .* cos_u );
    ddr_step =  2 .* dr_step .^2 ./ r_step - r_step .^2 .* ( 2 .* a2 - 2 .* a4 .* sin_u - ( a4 .* UIntervalStep + a3 ) .* cos_u + 2 .* a6 .* cos_u - ( a6 .* UIntervalStep + a5 ) .* sin_u );
    
    G_step = G(UIntervalStep) ;
    dG_step = G1 + G2 .* cos(UIntervalStep - u_0 + G3);
    
    F_r =  mu ./ r_step .^ 2 - G_step .^ 2 ./ r_step .^ 3 + (ddr_step .* G_step ./ r_step .^ 2 + dr_step .* dG_step ./ r_step .^ 2 - 2 .* dr_step .^ 2 .* dG_step ./ r_step .^ 3) .* G_step ./ r_step .^ 2 ;
    
    F_u = dG_step .* ( G_step ./ r_step .^ 3 ) ;
        

    uprime = 1./ (G_step./r_step.^2);
    TIME = TIME + (du/3)*(uprime(1) + 2*sum(uprime(3:2:end-1)) + ...
        + 4*sum(uprime(2:2:end-1)) + uprime(end));
    
    dvec = abs(sqrt(F_r.^2 + F_u.^2).* uprime);
    
    DV = DV + (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
        + 4*sum(dvec(2:2:end-1)) + dvec(end));
    

end



for ii = 1:ll
    
    if i_or_h(ii) == 0

        dU = U(ii,2) - U(ii,1);
        
        NumberInclinationSteps = dU ./ TotalUInclinationChange;  
        
        DeltaInclinationInThisStep = DeltaInclinationToFinalValue .* NumberInclinationSteps ;
        
        FinalIncStep = InitialIncStep + DeltaInclinationInThisStep;
        
        n_steps_i_h = ceil((U(ii,2)-U(ii,1))/du_i_h + 1) ;
        if n_steps_i_h < 5
            n_steps_i_h = 5 ;
        end
        UIntervalStep = linspace(U(ii,1),U(ii,2),n_steps_i_h)';
        
        sin_u = sin(UIntervalStep);
        cos_u = cos(UIntervalStep);
        
        r_step = r(UIntervalStep) ;
        dr_step = - r_step .^2 .* ( a1 + 2 .* a2 .* UIntervalStep + a4 .* cos_u - ( a4 .* UIntervalStep + a3 ) .* sin_u + a6 .* sin_u + ( a6 .* UIntervalStep + a5 ) .* cos_u );
        ddr_step =  2 .* dr_step .^2 ./ r_step - r_step .^2 .* ( 2 .* a2 - 2 .* a4 .* sin_u - ( a4 .* UIntervalStep + a3 ) .* cos_u + 2 .* a6 .* cos_u - ( a6 .* UIntervalStep + a5 ) .* sin_u );
        
        G_step = G(UIntervalStep) ;
        dG_step = G1 + G2 .* cos(UIntervalStep - u_0 + G3);
               
        du = UIntervalStep(2) - UIntervalStep(1) ;
        
        % Define boundary conditions
        u0 = UIntervalStep(1);
        i0 = InitialIncStep ;
        i1 = 0;
        i_aux = [dU^3 dU^2 ; 3*dU^2 2*dU]\[FinalIncStep-InitialIncStep ; 0 ] ;
        i2 = i_aux(2) ; 
        i3 = i_aux(1) ;
        
        DeltaU = (UIntervalStep - u0);
        i_step = i0 + i1 .* DeltaU + i2 .* DeltaU.^2 + i3 .* DeltaU .^ 3 ;
        di_step = i1 + 2.* i2 .* DeltaU + 3 .* i3 .* DeltaU .^ 2 ;
        ddi_step = 2.* i2 + 6 .* i3 .* DeltaU ;
        cot_i = cot(i_step);
        
        
        F_n =  G_step .^ 2 ./ r_step .^ 3 .* di_step ./ (cos_u + cot_i .* sin_u .* di_step);
        dF_n = ddi_step .* G_step .^ 2 ./ r_step .^ 3 ./ (cos_u + cot_i .* sin_u .* di_step) + 0.2e1 .* di_step .* G_step ./ r_step .^ 3 ./ (cos_u + cot_i .* sin_u .* di_step) .* dG_step - 0.3e1 .* di_step .* G_step .^ 2 ./ r_step .^ 4 ./ (cos_u + cot_i .* sin_u .* di_step) .* dr_step - di_step .* G_step .^ 2 ./ r_step .^ 3 ./ (cos_u + cot_i .* sin_u .* di_step) .^ 2 .* (-sin_u + di_step .^ 2 .* (-0.1e1 - cot_i .^ 2) .* sin_u + cot_i .* cos_u .* di_step + cot_i .* sin_u .* ddi_step);
        dF_n(1) = 0;
        
        F_r = mu ./ r_step .^ 2 - G_step .^ 2 ./ r_step .^ 3 + ((dG_step ./ (r_step .^ 2) - (2 .* G_step ./ r_step .^ 3 .* dr_step) - dr_step .* cot(i_step) .* sin_u ./ G_step .* F_n - r_step .* di_step .* (-0.1e1 - cot(i_step) .^ 2) .* sin_u ./ G_step .* F_n - r_step .* cot(i_step) .* cos_u ./ G_step .* F_n + r_step .* cot(i_step) .* sin_u ./ (G_step .^ 2) .* F_n .* dG_step - r_step .* cot(i_step) .* sin_u ./ G_step .* dF_n) .* dr_step + ((G_step ./ r_step .^ 2) - r_step .* cot(i_step) .* sin_u ./ G_step .* F_n) .* ddr_step) .* ((G_step ./ r_step .^ 2) - r_step .* cot(i_step) .* sin_u ./ G_step .* F_n);
        F_u = dG_step .* ( G_step ./ r_step .^ 3 - cot_i .* sin_u ./ G_step .* F_n );
        
        Hdot = r_step ./ G_step .* sin_u ./ sin(i_step) .* F_n ./ ( G_step ./ r_step .^2 - r_step .*cot_i .* sin_u ./ G_step .* F_n) ;
        FinalHStep = InitialHStep + (du/3)*(Hdot(1) + 2*sum(Hdot(3:2:end-1)) + ...
                  + 4*sum(Hdot(2:2:end-1)) + Hdot(end));

        uprime =  1./(G_step./r_step.^2-r_step.*cot_i.*sin_u.*F_n./G_step);
        dvec = abs(sqrt(F_r.^2 + F_u.^2 + F_n.^2).* uprime);
        
        DV = DV + (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
                  + 4*sum(dvec(2:2:end-1)) + dvec(end));
        TIME = TIME + (du/3)*(uprime(1) + 2*sum(uprime(3:2:end-1)) + ...
        + 4*sum(uprime(2:2:end-1)) + uprime(end));
        
        if U(ii,2) < u_f-1e-6

            if ii == ll
                
                n_steps = ceil((u_f - (U(ii,2)+du_others))/du_others +1) ;
                if n_steps < 5
                    n_steps = 5 ; 
                end
                UIntervalStep = linspace(U(ii,2),u_f,n_steps+1)';
                UIntervalStep = UIntervalStep(2:end);
                du = UIntervalStep(2)-UIntervalStep(1) ;
            else

                n_steps = ceil( ( (U(ii+1,1)-du_others) - (U(ii,2)+du_others) )/du_others +1 ) ;
                if n_steps < 5
                    n_steps = 5 ; 
                end
                UIntervalStep = linspace(U(ii,2),U(ii+1,1),n_steps+2)';
                UIntervalStep = UIntervalStep(2:end-1);
                du = UIntervalStep(2)-UIntervalStep(1) ;
            end
            
            sin_u = sin(UIntervalStep);
            cos_u = cos(UIntervalStep);
            
            r_step = r(UIntervalStep) ;
            dr_step = - r_step .^2 .* ( a1 + 2 .* a2 .* UIntervalStep + a4 .* cos_u - ( a4 .* UIntervalStep + a3 ) .* sin_u + a6 .* sin_u + ( a6 .* UIntervalStep + a5 ) .* cos_u );
            ddr_step =  2 .* dr_step .^2 ./ r_step - r_step .^2 .* ( 2 .* a2 - 2 .* a4 .* sin_u - ( a4 .* UIntervalStep + a3 ) .* cos_u + 2 .* a6 .* cos_u - ( a6 .* UIntervalStep + a5 ) .* sin_u );
            
            G_step = G(UIntervalStep) ;
            dG_step = G1 + G2 .* cos(UIntervalStep - u_0 + G3);
            
            F_r =  mu ./ r_step .^ 2 - G_step .^ 2 ./ r_step .^ 3 + (ddr_step .* G_step ./ r_step .^ 2 + dr_step .* dG_step ./ r_step .^ 2 - 2 .* dr_step .^ 2 .* G_step ./ r_step .^ 3) .* G_step ./ r_step .^ 2 ;
            
            F_u = dG_step .* ( G_step ./ r_step .^ 3 ) ;
            
            uprime =  1./(G_step./r_step.^2);
            dvec = abs(sqrt(F_r.^2 + F_u.^2).*uprime);
            
            DV = DV + (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
                  + 4*sum(dvec(2:2:end-1)) + dvec(end));
              
            TIME = TIME + (du/3)*(uprime(1) + 2*sum(uprime(3:2:end-1)) + ...
                  + 4*sum(uprime(2:2:end-1)) + uprime(end));
                      
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
        
        n_steps_i_h = ceil((U(ii,2)-U(ii,1))/du_i_h + 1) ;
        if n_steps_i_h < 5 % 4-degree polyfit
            n_steps_i_h = 5 ;
        end
        UIntervalStep = linspace(U(ii,1),U(ii,2),n_steps_i_h)';
        
        sin_u = sin(UIntervalStep);
        cos_u = cos(UIntervalStep);
        
        r_step = 1./( a0 + a1.* UIntervalStep + a2.* UIntervalStep .^2 + ( a3 + a4.* UIntervalStep ).*cos(UIntervalStep) + ( a5 + a6.*UIntervalStep ).*sin_u ) ;
        dr_step = - r_step .^2 .* ( a1 + 2 .* a2 .* UIntervalStep + a4 .* cos_u - ( a4 .* UIntervalStep + a3 ) .* sin_u + a6 .* sin_u + ( a6 .* UIntervalStep + a5 ) .* cos_u );
        ddr_step =  2 .* dr_step .^2 ./ r_step - r_step .^2 .* ( 2 .* a2 - 2 .* a4 .* sin_u - ( a4 .* UIntervalStep + a3 ) .* cos_u + 2 .* a6 .* cos_u - ( a6 .* UIntervalStep + a5 ) .* sin_u );
        
        G_step = G0 + G1.*(UIntervalStep-u_0) + G2.*sin( UIntervalStep-u_0 + G3 );
        dG_step = G1 + G2 .* cos(UIntervalStep - u_0 + G3);
        
        du = UIntervalStep(2) - UIntervalStep(1);
                    
        % Define boundary conditions
        
        u0 = UIntervalStep(1);
        h0 = InitialHStep ;
        h1 = 0;
        h_aux = [dU^3 dU^2 ; 3*dU^2 2*dU]\[FinalHStep-InitialHStep ; 0] ;
        h2 = h_aux(2) ;
        h3 = h_aux(1) ;
        
        DeltaU = (UIntervalStep - u0);
        H_Step = h0 + h1 .* DeltaU + h2 .* DeltaU .^ 2 + h3 .* DeltaU .^ 3 ;
        dH_step = h1 + 2 .* h2 .* DeltaU + 3 .* h3 .* DeltaU .^ 2 ;
        ddH_step = 2 .* h2 + 6 .* h3 .* DeltaU ;
    
        F_n_Step = @(u) (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_0) + G2 .* sin(u - u_0 + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(InitialIncStep) + cot(InitialIncStep) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2));       
        iStep_old = InitialIncStep .* ones(length(UIntervalStep),1);
        
        %% PUOI PENSARE A UN ITERATIVE METHOD
        
%         Idot = r_step ./ G_step .* cos_u .* F_n ./ ( G_step ./ r_step .^ 2 - r_step .* cot(x) .*sin(u)/G(u) .* F_n(u) ) 
%         t1 = tic;
%         sol = ode45(@(u,x) DiffInclinationHEquationsWithControl(u,x,r,G,F_n_Step,0),UIntervalStep,InitialIncStep);      
%         [i_step] = deval(sol,UIntervalStep)';
%         tf1 = toc(t1)
%         

        F_n = dH_step .* G_step .^ 2 ./ r_step .^ 3 ./ (sin_u ./ sin(iStep_old) + cot(iStep_old) .* sin_u .* dH_step);
        di_step = r_step ./ G_step .* cos_u .* F_n ./ ( G_step ./ r_step .^2 - r_step .*cot(iStep_old) .* sin_u ./ G_step .* F_n) ;
        i_step = zeros(n_steps_i_h,1) ;
        i_step(1) = InitialIncStep ;
        for K = 2 : n_steps_i_h
            i_step(K) = i_step(K-1) + 0.5*(di_step(K-1)+di_step(K))*du ;
        end
        
        
        jj = 0 ;

        while max(abs(i_step-iStep_old))>1e-6 && jj < 5
            jj = jj +1 ;
            iStep_old = i_step;
            
%             ti1 = tic;
%             ppp = polyfit(UIntervalStep-u0,iStep_old,4);
%             %         INC = @(u) ppp(1).*(u-u0).^6+ppp(2).*(u-u0).^5+ppp(3).*(u-u0).^4+ppp(4).*(u-u0).^3+ppp(5).*(u-u0).^2+ppp(6).*(u-u0)+ppp(7);
%             INC = @(u) ppp(1).*(u-u0).^4+ppp(2).*(u-u0).^3+ppp(3).*(u-u0).^2+ppp(4).*(u-u0)+ppp(5);
%             %   Interpolo con degree 4
%             F_n_Step = @(u) (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2) .* (G0 + G1 .* (u - u_0) + G2 .* sin(u - u_0 + G3)) .^ 2 .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 3 ./ (sin(u) ./ sin(INC(u)) + cot(INC(u)) .* sin(u) .* (h1 + 2 .* h2 .* (u - u0) + 3 .* h3 .* (u - u0) .^ 2));
%             sol = ode45(@(u,x) DiffInclinationHEquationsWithControl(u,x,r,G,F_n_Step,0),UIntervalStep,InitialIncStep);
%             [i_step] = deval(sol,UIntervalStep)';
%             tf1 = toc(ti1)
            
            F_n = dH_step .* G_step .^ 2 ./ r_step .^ 3 ./ (sin_u ./ sin(iStep_old) + cot(iStep_old) .* sin_u .* dH_step);
            di_step = r_step ./ G_step .* cos_u .* F_n ./ ( G_step ./ r_step .^2 - r_step .*cot(iStep_old) .* sin_u ./ G_step .* F_n) ;
            
            i_step = zeros(n_steps_i_h,1) ;
            i_step(1) = InitialIncStep ;
            for K = 2 : n_steps_i_h
            i_step(K) = i_step(K-1) + 0.5*(di_step(K-1)+di_step(K))*du ;
            end
            
        end

        if jj == 5
            warning('Non interpolato bene')
        end
        
        FinalIncStep = i_step(end) ;
        sin_i = sin(i_step) ;
        cot_i = cot(i_step) ;
        
%         if jj~=0
%             di_step = 4.*ppp(1).*(UIntervalStep-u0).^3+3.*ppp(2).*(UIntervalStep-u0).^2+2.*ppp(3).*(UIntervalStep-u0)+ppp(4) ;
%           
%         else
%             
%             ppp = polyfit(UIntervalStep-u0,i_step,1);
%             di_step = ppp(1);
%         end
        
%         F_n =  F_n_Step(UIntervalStep);
        dF_n = ddH_step .* G_step .^ 2 ./ r_step .^ 3 ./ (sin_u ./ sin_i + cot_i .* sin_u .* dH_step) + 0.2e1 .* dH_step .* G_step ./ r_step .^ 3 ./ (sin_u ./ sin_i + cot_i .* sin_u .* dH_step) .* dG_step - 0.3e1 .* dH_step .* G_step .^ 2 ./ r_step .^ 4 ./ (sin_u ./ sin_i + cot_i .* sin_u .* dH_step) .* dr_step - dH_step .* G_step .^ 2 ./ r_step .^ 3 ./ (sin_u ./ sin_i + cot_i .* sin_u .* dH_step) .^ 2 .* (cos_u ./ sin_i - sin_u ./ sin_i .^ 2 .* di_step .* cos(i_step) + di_step .* (-0.1e1 - cot_i .^ 2) .* sin_u .* dH_step + cot_i .* cos_u .* dH_step + cot_i .* sin_u .* ddH_step);
        dF_n(1) = 0;

%         i_step = INC(UIntervalStep) ;
      
        F_r = mu ./ r_step .^ 2 - G_step .^ 2 ./ r_step .^ 3 + (ddr_step .* (G_step ./ r_step .^ 2 - r_step .* cot_i .* sin_u ./ G_step .* F_n ) + dr_step .* (dG_step ./ r_step .^ 2 - 0.2e1 .* G_step ./ r_step .^ 3 .* dr_step - dr_step .* cot_i .* sin_u ./ G_step .* F_n - r_step .* di_step .* (-0.1e1 - cot_i .^ 2) .* sin_u ./ G_step .* F_n - r_step .* cot_i .* cos_u ./ G_step .* F_n + r_step .* cot_i .* sin_u ./ G_step .^ 2 .* F_n .* dG_step - r_step .* cot_i .* sin_u ./ G_step .* dF_n)) .* (G_step ./ r_step .^ 2 - r_step .* cot_i .* sin_u ./ G_step .* F_n);
        
        F_u = dG_step .* ( G_step ./ r_step .^ 3 - cot_i .* sin_u ./ G_step .* F_n );
       
        uprime =  1./(G_step./r_step.^2-r_step.*cot_i.*sin_u.*F_n./G_step);
        dvec = abs(sqrt(F_r.^2 + F_u.^2 + F_n.^2).* uprime);
        
        DV = DV + (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
                  + 4*sum(dvec(2:2:end-1)) + dvec(end));
              
        TIME = TIME + (du/3)*(uprime(1) + 2*sum(uprime(3:2:end-1)) + ...
               + 4*sum(uprime(2:2:end-1)) + uprime(end));
                
        if U(ii,2) < u_f-1e-6

            if ii == ll
                
                n_steps = ceil((u_f - (U(ii,2)+du_others))/du_others +1) ;
                if n_steps < 5
                    n_steps = 5 ; 
                end
                UIntervalStep = linspace(U(ii,2),u_f,n_steps+1)';
                UIntervalStep = UIntervalStep(2:end);
                du = UIntervalStep(2)-UIntervalStep(1) ;
                
            else
                
                n_steps = ceil( ( (U(ii+1,1)-du_others) - (U(ii,2)+du_others) )/du_others +1 ) ;
                if n_steps < 5
                    n_steps = 5 ; 
                end
                UIntervalStep = linspace(U(ii,2),U(ii+1,1),n_steps+2)';
                UIntervalStep = UIntervalStep(2:end-1);
                du = UIntervalStep(2)-UIntervalStep(1) ;
            end
            
            sin_u = sin(UIntervalStep);
            cos_u = cos(UIntervalStep);
            
            r_step = 1./( a0 + a1.* UIntervalStep + a2.* UIntervalStep .^2 + ( a3 + a4.* UIntervalStep ).*cos(UIntervalStep) + ( a5 + a6.*UIntervalStep ).*sin(UIntervalStep) ) ;
            dr_step = - r_step .^2 .* ( a1 + 2 .* a2 .* UIntervalStep + a4 .* cos_u - ( a4 .* UIntervalStep + a3 ) .* sin_u + a6 .* sin_u + ( a6 .* UIntervalStep + a5 ) .* cos_u );
            ddr_step =  2 .* dr_step .^2 ./ r_step - r_step .^2 .* ( 2 .* a2 - 2 .* a4 .* sin_u - ( a4 .* UIntervalStep + a3 ) .* cos_u + 2 .* a6 .* cos_u - ( a6 .* UIntervalStep + a5 ) .* sin_u );
            
            G_step = G0 + G1.*(UIntervalStep-u_0) + G2.*sin( UIntervalStep-u_0 + G3 );
            dG_step = G1 + G2 .* cos(UIntervalStep - u_0 + G3);
            
            F_r =  mu ./ r_step .^ 2 - G_step .^ 2 ./ r_step .^ 3 + (ddr_step .* G_step ./ r_step .^ 2 + dr_step .* dG_step ./ r_step .^ 2 - 2 .* dr_step .^ 2 .* G_step ./ r_step .^ 3) .* G_step ./ r_step .^ 2 ;
            
            F_u = dG_step .* ( G_step ./ r_step .^ 3 ) ;
            
%             History = [History ; UIntervalStep ] ;
            
%             d_step = [d_step ; abs(sqrt(F_r.^2 + F_u.^2)./(G_step./r_step.^2))] ;
            uprime =  1./(G_step./r_step.^2);
            dvec = abs(sqrt(F_r.^2 + F_u.^2).* uprime);
            
            DV = DV + (du/3)*(dvec(1) + 2*sum(dvec(3:2:end-1)) + ...
                + 4*sum(dvec(2:2:end-1)) + dvec(end));
            TIME = TIME + (du/3)*(uprime(1) + 2*sum(uprime(3:2:end-1)) + ...
                  + 4*sum(uprime(2:2:end-1)) + uprime(end));
%             FFF = [FFF ; sqrt(F_r.^2 + F_u.^2)];
                
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
    
