% function  [I, dI] = Integrals_Drag3_2(L0, L, P10, P20, a0, drag, R, e_flag)
function  [I, dI] = Integrals_Drag3_2(L0, L, P10, P20, Q10, Q20, a0, drag, R, J2, e_flag, Earth_flat)


% Calculates analytical integrals for the drag acceleration
%
% INPUTS:
% L0: longitude of the initial state.
% L: longitude of the current state.
% P10, P20: values at x(1) for Equinoctial elements P1 and P2.
% a0: semimajor axis
% drag: structure with the input for the drag
% e_flag:
% drag.num_an: flag. 1 to compute integrals analytically, 0 to compute them
% numerically
%
% OUTPUTS:
% I: 6-row vector with the integrals between x(1) and x(2).


% Marilena Di Carlo, 2015



global checking_crossed_intervals checking_altitudes
checking_crossed_intervals = [];
checking_altitudes = [];
%% Constants - dargliele come input!



% % Adimensional Earth Radius
% R = 1;



if any(abs(L-L0)>1e-6)
    
    [a,b] = size(L);
    
    if a > b
        L=L.';
    end
    
    e0 = sqrt(P10^2 + P20^2);
    
    %% Useful precomputation - copied by Integrals_tang.m
    
    c_L0 = cos(L0);
    s_L0 = sin(L0);
    
    % e*cos(theta0) where L = (theta + Omega + omega)
    e0_c_th0 = P10*s_L0+P20*c_L0;
    
    % e*sin(theta0) where L = (theta + Omega + omega)
    e0_s_th0 = P20*s_L0-P10*c_L0;
    
    % theta0
    th0 = atan2(e0_s_th0,e0_c_th0);
    
    % -----------------------------------------------------------------
    % Theta computed as Federico did - I got singularity so I replaced
    % it with another method, below
    % -----------------------------------------------------------------
    
    %         % e*cos(theta) where L = (theta + Omega + omega)
    %         e0_c_th = P10*s_L+P20*c_L;
    %
    %         % e*sin(theta) where L = (theta + Omega + omega)
    %         e0_s_th = P20*s_L-P10*c_L;
    %
    %         % theta
    %         th = atan2(e0_s_th,e0_c_th);
    
    % -----------------------------------------------------------------
    % new method to compute theta
    % -----------------------------------------------------------------
    
    % Omega + omega
    Omom = atan2(P10, P20);
    
    % Omega
    Omega = atan2(Q10, Q20);
    
    % omega
    omega = mod(Omom - Omega, 2*pi);
    
    % inclination
    incl = 2 * atan(sqrt(Q10^2 + Q20^2));
    
    % theta
    th = L - Omom;
    
    th0 = L0 - Omom;
    
    
    %% Find in which height range for the interpolation of the density the orbit takes place
    
    % Semiparameter
    p = a0 * (1-e0^2);
    
    % h_crossing defines the height limits for each region in which the
    % coefficients of the Chebyshev expansions are the same ones
    h_crossing = drag.h_crossing;    
    
    % Apogee and perigee height [DU]
    h_perigee = a0 * (1-e0) - R;
    h_apogee  = a0 * (1+e0) - R;
    
    % Identify region of the Chebyshev fitting which are crossed by the orbit
    [~,index1] = find(h_crossing < h_perigee);
    [~,index2] = find(h_crossing > h_apogee);
    
    % Define vector h_crossing for the current orbit, h_crossing_orbit
    % 0 element vector means that the orbit lies entirely inside one region
    % 1 element vector means that the orbit lies in two regions
    % 2 elements vector means that the orbit lies in three regions
    % .... and so on
    if ~isempty(index1) && ~isempty(index2)
        h_crossing_orbit = h_crossing(:, index1(end) + 1 : index2 - 1);
    else
        warning('Orbit outside boundaries for the interpolation of the density')
        disp(strcat('Apogee altitude is [km]: ',num2str(h_apogee*6378.136)))
        disp(strcat('Perigee altitude is [km]: ',num2str(h_perigee*6378.136)))
        keyboard
    end
    
    
    % Theta corresponding to the crossing conditions from one region to
    % another
    theta_crossing_orbit = acos((1/e0) .*( (p./(h_crossing_orbit + R)) -1));
    
    % Add also 2*pi - theta_crossing_orbit to consider the entire theta
    % range for the orbit, from 0 to 2pi
    theta_crossing_orbit0 = [theta_crossing_orbit 2*pi-theta_crossing_orbit];
    theta_crossing_orbit = theta_crossing_orbit0;

%     keyboard
    if ~isempty(theta_crossing_orbit)
        n = 1;
        while theta_crossing_orbit(end) < th(end)
            theta_crossing_orbit = [theta_crossing_orbit (theta_crossing_orbit0 + 2 *n*pi)];
            n = n+1;
        end
        
        n = 1;
        while min(theta_crossing_orbit) > th0
            theta_crossing_orbit = [theta_crossing_orbit (theta_crossing_orbit0 - 2 *n*pi)];
            n = n+1;
        end
    end
 
    theta_crossing_orbit = sort(theta_crossing_orbit);

    
    
    
    %% Integrals computation
    
    % Initial point height [DU]
    h0 = a0 * (1-e0^2) ./ (1+e0*cos(th0)) - R;
    
    % In which interval for the definition of the Chebyshev expansion does
    % the initial point falls?
    for k = 1 : length(h_crossing)-1
        
        if h0 < h_crossing(k+1) && h0 > h_crossing(k)
            interval0 = k;
            break
        end
        
%         % Correggere! 
%         if h0 == h_crossing(k)
%             interval0 = k;
%         end
    end
    
%     if ~exist('interval0')
%         keyboard
%     end
    
    % Height corresponding to each angle th
    h = a0 * (1-e0^2) ./ (1+e0*cos(th)) - R;
    
    % In which interval for the definition of the Chebyshev expansion do
    % these point falls?
%     interval = zeros(1,length(th));
    for i = 1 : length(th)
        
        for k = 1 : length(h_crossing) - 1
            
            if h(i) < h_crossing(k+1) && h(i) > h_crossing(k)
                interval(1,i) = k;
                break
            end
        end
        
    end
    
    % Current integration may not consider entire orbit. Consider only
    % value of theta_crossing greater than th0 and lower than th(end)
    theta_crossing_interval = theta_crossing_orbit(theta_crossing_orbit>th0);
    theta_crossing_interval = theta_crossing_interval(theta_crossing_interval<th(end));
    
    
    % The variable all_theta is very important to define the region where
    % the Chebychev polynomials are defined:
    % 1st row: true longitude value theta (both the one defined by th and
    %          the one when a crossin of a region occurs)
    % 2nd row: number identifying the region for the definition of the
    %          polynomials (0 if the angle is a crossing angle)
    % 3rd row: height corresponding to each angle. NEcessary when there are
    %          multiple crossing theta close to each other for which the region
    %          needs to be defined
    % 4th row: number identifying the region for the integration of the
    %          Chebyshev polynomials FOR THE CROSSING ANGLES ONLY
    
    % ---------------------------------------------------------------------
    % Non sono sicura che la seguente prima condizione if serva davvero...
    % Commentata per ora, nel caso decommentare
    % ---------------------------------------------------------------------
%     if th(1,end)>2*pi
%         all_theta(1,:) = [th0       theta_crossing_interval          2*pi-1e-6  2*pi      th];
%         all_theta(2,:) = [interval0 zeros(1,length(theta_crossing_interval)+2) interval];
%     else
%         all_theta(1,:) = th0;
%         all_theta(1,:) = [all_theta(1,:) th];
         all_theta(1,:) = [th0       theta_crossing_interval                  th];
         all_theta(2,:) = [interval0 zeros(1,length(theta_crossing_interval)) interval];
%     end
    

    all_theta(3,:) = a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,:))) - R;
    all_theta(4,:) = zeros(1,length(all_theta(3,:)));

    all_theta = sortrows(all_theta',1)';


    % Initialization 
    % Total integrals, given by the sum of IDrag1, IDrag2 etc multiplied by
    % appropriate coefficients
    I_drag1 = zeros(1,length(th));
    I_drag2 = zeros(1,length(th));
    I_drag4 = zeros(1,length(th));
    I_drag5 = zeros(1,length(th));
    I_drag6 = zeros(1,length(th));
    
    % Temporary value of the integral before saving of the solution
    I_drag1_tmp = 0;
    I_drag2_tmp = 0;
    I_drag4_tmp = 0;
    I_drag5_tmp = 0;
    I_drag6_tmp = 0;
    
    % Initial integration angle
    theta_in = th0;
    
    % Initial index to save solutions
    save_index = 1;
    
    % Computed analytical integrals - single small integrals - see
    % references
%     keyboard
    [IDrag1, IDrag2, IDrag4, IDrag5, IDrag6] = Integrals_Drag_Chebyshev3_2(e0, theta_in, all_theta(1,2:end));

    % Now we need to multiply the integrals IDrag1, IDrag2 etc by
    % appropriate coefficients that depend upon the Chebyshev coefficients,
    % and therefore upon the range of height in which we are for every
    % angle theta
    % The following is to determine where we are. 
    for  i = 2 : length(all_theta(1,:))

        % Final integration angle
        theta_f = all_theta(1,i);
        
        
        
        % Which coefficients to use in the integals?
        
        % ---------------------------------------------------------------
%       % Questo non dovrebbe mai succedere credo.... eliminare questa
%       % condizione!!!!!!!!!
        % ---------------------------------------------------------------
%         if all_theta(1,i) == 2*pi-0.0001
%             
%             all_theta(2,i) = max(all_theta(2,i-1),all_theta(4,i-1));
%             coeff_current = drag.coeff(:,all_theta(2,i));
          
          % ---------------------------------------------------------------
%         % Questo non dovrebbe mai succedere credo.... eliminare questa
%         % condizione!!!!!!!!!
          % ---------------------------------------------------------------
%         elseif all_theta(1,i) == 0   && i~=1
% 
%             all_theta(2,i) = max(all_theta(2,i-1),all_theta(4,i-1));
%             coeff_current = drag.coeff(:,all_theta(2,i));
            
        % If a number if defined in the second row of all_theta, that's
        % the region to use!
        if all_theta(2,i) ~= 0
%             keyboard
            coeff_current =  drag.coeff(:,all_theta(2,i));
            checking_crossed_intervals = [checking_crossed_intervals; all_theta(2,i)];
            checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
            
        % If the current angle is a crossing angle but the previous
        % one was not, propagate from the previous to the current
        % using the region identified by the previous
        elseif all_theta(2,i) == 0  && all_theta(2,i-1) ~=0
            
            all_theta(4,i) = all_theta(2,i-1);
            coeff_current =  drag.coeff(:,all_theta(2,i-1));
            checking_crossed_intervals = [checking_crossed_intervals; all_theta(2,i-1)];
            checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
%             keyboard
        % If the current angle is a crossing angle and also the
        % previous angle was a crossing angle, than everything
        % depends on what happens to the previous previous
        % angle....
        elseif all_theta(2,i) == 0 && all_theta(2,i-1) == 0
%             keyboard
            % If the previous previous angle was not a crossing angle
            if all_theta(2,i-2) ~= 0
                
                % If height is increasing, use the region of the
                % previous previous angle plus 1
                if all_theta(3,i) > all_theta(3,i-1)
                    all_theta(4,i) = all_theta(2,i-2)+1;
                    coeff_current = drag.coeff(:,all_theta(2,i-2)+1);
                    checking_crossed_intervals = [checking_crossed_intervals; all_theta(2,i-2)+1];
                    checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
                    
                % If height is decreasing use region of the
                % previous previous angle minus 1
                elseif all_theta(3,i) < all_theta(3,i-1)
                    
                    all_theta(4,i) = all_theta(2,i-2)-1;
                    coeff_current = drag.coeff(:,all_theta(2,i-2)-1);
                    checking_crossed_intervals = [checking_crossed_intervals; all_theta(2,i-2)-1];
                    checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
                end
                
            % If the previous previous angle was also a crossing
            % angle, use the 4th row of the variable all_theta to
            % understand what to do
            else
%                 keyboard
                if all_theta(3,i) > all_theta(3,i-1)
                    
                    all_theta(4,i) = all_theta(4,i-1)+1;
                    coeff_current = drag.coeff(:,all_theta(4,i-1)+1);
                    checking_crossed_intervals = [checking_crossed_intervals; all_theta(4,i-1)+1];
                    checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
                    
                elseif all_theta(3,i) < all_theta(3,i-1)

                    all_theta(4,i) = all_theta(4,i-1)-1;
                    coeff_current = drag.coeff(:,all_theta(4,i-1)-1);
                    checking_crossed_intervals = [checking_crossed_intervals; all_theta(4,i-1)-1];
                    checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];
                    
                elseif all_theta(3,i) == all_theta(3,i-1)
                    
                    all_theta(4,i) = all_theta(4,i-1);
                    coeff_current = drag.coeff(:,all_theta(4,i-1));
                    checking_crossed_intervals = [checking_crossed_intervals; all_theta(4,i-1)];
                    checking_altitudes = [checking_altitudes; (a0 * (1-e0^2) ./ (1+e0*cos(all_theta(1,i))) - R)*6378.136];

                end
            end

        end
        

        if drag.num_an == 1
            
            % +============================================================
            % Analytical Integrals
            % -------------------------------------------------------------
            % Coefficients of the Chebyshev expansion are saved into
            % coeff_currents as:
            % coeff_current = [a8; a7; a6; a5; a4; a3; a2; a1; a0];
            % They are flipped here for ease of use
            coeff = flipud(coeff_current);
            
            
            % Rendere k automatizzabile in base all'ordine
            % dell'espansione
            if length(coeff) == 8
                k = [( coeff(1)     - coeff(2)*R  +   coeff(3)*R^2 -   coeff(4)*R^3   +   coeff(5)*R^4     -    coeff(6)*R^5     +    coeff(7)*R^6     -    coeff(8)*R^7     ); ...
                    (                coeff(2)*p  - 2*coeff(3)*p*R + 3*coeff(4)*p*R^2 - 4*coeff(5)*p*R^3   +  5*coeff(6)*p*R^4   -  6*coeff(7)*p*R^5   +  7*coeff(8)*p*R^6   );  ...
                    (                                coeff(3)*p^2 - 3*coeff(4)*p^2*R + 6*coeff(5)*p^2*R^2 - 10*coeff(6)*p^2*R^3 + 15*coeff(7)*p^2*R^4 - 21*coeff(8)*p^2*R^5 ); ...
                    (                                                 coeff(4)*p^3   - 4*coeff(5)*p^3*R   + 10*coeff(6)*p^3*R^2 - 20*coeff(7)*p^3*R^3 + 35*coeff(8)*p^3*R^4 ); ...
                    (                                                                    coeff(5)*p^4     -  5*coeff(6)*p^4*R   + 15*coeff(7)*p^4*R^2 - 35*coeff(8)*p^4*R^3 ); ...
                    (                                                                                          coeff(6)*p^5     -  6*coeff(7)*p^5*R   + 21*coeff(8)*p^5*R^2 ); ...
                    (                                                                                                                coeff(7)*p^6     -  7*coeff(8)*p^6*R   ); ...
                    (                                                                                                                                      coeff(8)*p^7     )] ;
                
            elseif length(coeff) == 6
                k = [( coeff(1)     - coeff(2)*R  +   coeff(3)*R^2 -   coeff(4)*R^3   +   coeff(5)*R^4     -    coeff(6)*R^5          ); ...
                    (                coeff(2)*p  - 2*coeff(3)*p*R + 3*coeff(4)*p*R^2 - 4*coeff(5)*p*R^3   +  5*coeff(6)*p*R^4    );  ...
                    (                                coeff(3)*p^2 - 3*coeff(4)*p^2*R + 6*coeff(5)*p^2*R^2 - 10*coeff(6)*p^2*R^3 ); ...
                    (                                                 coeff(4)*p^3   - 4*coeff(5)*p^3*R   + 10*coeff(6)*p^3*R^2 ); ...
                    (                                                                    coeff(5)*p^4     -  5*coeff(6)*p^4*R   ); ...
                    (                                                                                          coeff(6)*p^5    )] ;
                
                
            elseif length(coeff) == 5
                k = [( coeff(1)     - coeff(2)*R  +   coeff(3)*R^2 -   coeff(4)*R^3   +   coeff(5)*R^4       ); ...
                    (                coeff(2)*p  - 2*coeff(3)*p*R + 3*coeff(4)*p*R^2 - 4*coeff(5)*p*R^3   );  ...
                    (                                coeff(3)*p^2 - 3*coeff(4)*p^2*R + 6*coeff(5)*p^2*R^2  ); ...
                    (                                                 coeff(4)*p^3   - 4*coeff(5)*p^3*R   ); ...
                    (                                                                    coeff(5)*p^4    )] ;
                

            end

            I_drag1_tmp = I_drag1_tmp + dot(k,IDrag1(1:length(coeff),i-1));
            I_drag2_tmp = I_drag2_tmp + dot(k,IDrag2(1:length(coeff),i-1));
            I_drag4_tmp = I_drag4_tmp + dot(k,IDrag4(1:length(coeff),i-1));
            I_drag5_tmp = I_drag5_tmp + dot(k,IDrag5(1:length(coeff),i-1));
            I_drag6_tmp = I_drag6_tmp + dot(k,IDrag6(1:length(coeff),i-1));

       else
           % +============================================================
            % Numerical Integrals
            % -------------------------------------------------------------
            
%             I_drag1_tmp =  I_drag1_tmp + integral(@(theta)dIdrag1_dL3(theta,a0,e0, coeff_current),theta_in,theta_f);
%             I_drag2_tmp =  I_drag2_tmp + integral(@(theta)dIdrag2_dL3(theta,a0,e0, coeff_current),theta_in,theta_f);
%             I_drag4_tmp =  I_drag4_tmp + integral(@(theta)dIdrag4_dL3(theta,a0,e0, coeff_current),theta_in,theta_f);
%             I_drag5_tmp =  I_drag5_tmp + integral(@(theta)dIdrag5_dL3(theta,a0,e0, coeff_current),theta_in,theta_f);
%             I_drag6_tmp =  I_drag6_tmp + integral(@(theta)dIdrag6_dL3(theta,a0,e0, coeff_current),theta_in,theta_f);
            
            I_drag1_tmp =  I_drag1_tmp + integral(@(theta)dIdrag1_dL3(theta,a0,e0, incl, omega, coeff_current, R, J2, Earth_flat),theta_in,theta_f);
            I_drag2_tmp =  I_drag2_tmp + integral(@(theta)dIdrag2_dL3(theta,a0,e0, incl, omega, coeff_current, R, J2, Earth_flat),theta_in,theta_f);
            I_drag4_tmp =  I_drag4_tmp + integral(@(theta)dIdrag4_dL3(theta,a0,e0, incl, omega, coeff_current, R, J2, Earth_flat),theta_in,theta_f);
            I_drag5_tmp =  I_drag5_tmp + integral(@(theta)dIdrag5_dL3(theta,a0,e0, incl, omega, coeff_current, R, J2, Earth_flat),theta_in,theta_f);
            I_drag6_tmp =  I_drag6_tmp + integral(@(theta)dIdrag6_dL3(theta,a0,e0, incl, omega, coeff_current, R, J2, Earth_flat),theta_in,theta_f);
            
        end
        
        % Save integral only if point is not a crossing
        if all_theta(2,i) ~= 0
            
            I_drag1(1,save_index) = I_drag1_tmp;
            I_drag2(1,save_index) = I_drag2_tmp;
            I_drag4(1,save_index) = I_drag4_tmp;
            I_drag5(1,save_index) = I_drag5_tmp;
            I_drag6(1,save_index) = I_drag6_tmp;
            
            save_index = save_index + 1;
            
        end

        theta_in = theta_f;
            
    end

    I_drag7 = I_drag4 - I_drag1;



    
    I = [I_drag1; I_drag2; 0*I_drag2; I_drag4; I_drag5; I_drag6; I_drag7];
    dI = [0; 0; 0; 0; 0;0 ;0];
    
    
else
    
    I = zeros(7,length(L));
    dI = zeros(7,1);
    
end


% if ~exist('I_drag1','var')
%     warning('I_drag1 not defined')
%     keyboard
% end
% if abs(I_drag1(1,end)) > 1e14
%     keyboard
% end

end




