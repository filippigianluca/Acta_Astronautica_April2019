function coefficients = Chebyshev_interpolation(f, n, a, b, flag_log)

% Function for the interpolation of a function with Chebyshev
% Input: f -> function handle
%        n -> order of the polynomial
%        a -> lower boundary of the interval
%        b -> upper boundary of the interval
%        flag_log -> 1 if working with the logarithm of the function

% Marilena Di Carlo, 2015


% Compute nodes for Chebyshev and coefficient
[xi, c] = chebpolfit((f),n,a,b,flag_log);

% Define point for the polynomial 
x = linspace(a,b,1001)';

% Compute polynomial
P = chebpolval(x, c, n, a, b);

% Real values of the function at x and at the Chebyshev nodes
real = zeros(1,length(x));
f_nodes = zeros(1,length(xi));

if flag_log == 1
    % Real values at x
    for i = 1 : length(x)
        real(i) = log(feval(f,x(i)));
    end
    % Real values at Chebyshev nodes
    for k = 1 : length(xi)
        f_nodes(k) = log(feval(f,xi(k)));
    end
else
    % Real values at x
    for i = 1 : length(x)
        real(i) = (feval(f,x(i)));
    end
    % Real values at Chebyshev nodes
    for k = 1 : length(xi)
        f_nodes(k) = (feval(f,xi(k)));
    end
end


if flag_log == 1
    figure
    subplot(4,1,1)
    plot(x,real,'b','LineWidth',2)
    hold on
    plot(x,P,'r','LineWidth',2)
    hold on
    plot(xi,f_nodes,'ok','MarkerFaceColor',[0 0 0])
    % xlabel('h [km]')
%     ylabel('log {\rho} [kg/m^3]')
    ylabel('log()')
    legend('Real','Chebyshev','Chebyshev nodes')
    subplot(4,1,2)
    plot(x,(real-P'))
    % xlabel('h [km]')
    ylabel('log(real) - log(Chebyshev)')
    subplot(4,1,3)
    plot(x,exp(real),'b','LineWidth',2)
    hold on
    plot(x,exp(P),'r','LineWidth',2)
    hold on
    plot(xi,exp(f_nodes),'ok','MarkerFaceColor',[0 0 0])
    % xlabel('h [km]')
    % ylabel('{\rho} [kg/m^3]')
    subplot(4,1,4)
    plot(x,exp(real)-exp(P'))
    % xlabel('h [km]')
    ylabel('real - Chebyshev')
else 
    figure
    subplot(2,1,1)
    plot(x,real,'b','LineWidth',2)
    hold on
    plot(x,P,'r','LineWidth',2)
    hold on
    plot(xi,f_nodes,'ok','MarkerFaceColor',[0 0 0])
%     xlabel('h [km]')
%     ylabel('Chebyshev {\rho} [kg/m^3]')
    legend('Real Chebyshev','Chebyshev Chebysev','Chebyshev nodes')
    subplot(2,1,2)
    plot(x,(real-P'))
    ylabel('real - Chebyshev')
%     xlabel('h [km]')
%     ylabel('Density difference [kg/m^3]')
end


P_final = symbolic_result(c, n, a, b);
coefficients = sym2poly(P_final);



% %% Verifica
% n = length(coefficients);
% 
% if flag_log == 1
%     
%     for i = 1 : length(x)
%         rho_finale(i) = 0;
%         
%         for k = 1 : n
%             
%             rho_finale(i) = rho_finale(i) + coefficients(n-k+1) * x(i).^(k-1);
%             
%             
%         end
%         rho_reale(i) = log(exponential_atm_model(x(i)));
%     end
%     
% else
%     
%         for i = 1 : length(x)
%         rho_finale(i) = 0;
%         
%         for k = 1 : n
%             
%             rho_finale(i) = rho_finale(i) + coefficients(n-k+1) * x(i).^(k-1);
%             
%             
%         end
%         rho_reale(i) = log(exponential_atm_model(x(i)));
%     end
%     
% end
% 
% 
% figure
% plot(x,(rho_finale),'r')
% hold on
% plot(x,rho_reale,'b')
% xlabel('h [km]')
% ylabel('log {\rho} [kg/m^3]')
% legend('Chebyshev \rho','Real \rho')
% 
% 
% figure
% plot(x,exp(rho_finale),'r')
% hold on
% plot(x,exp(rho_reale),'b')
% xlabel('h [km]')
% ylabel('{\rho} [kg/m^3]')
% legend('Chebyshev \rho','Real \rho')
