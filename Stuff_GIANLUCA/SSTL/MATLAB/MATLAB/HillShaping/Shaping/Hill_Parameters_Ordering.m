function [parameters] = Hill_Parameters_Ordering (x,l,shape_flag,param_flag)

% Hill_Parameters_Ordering
% Function that orders the BC parameters and the free-parameters in one
% vector, according to their appearence in the shaping functions
%
% INPUTs: 
% ...
% SCRIVIIIIIIIIIIIIIIIIIIIIII

if (shape_flag == 0) % Radial: 1/[a0+a1*u+a2*u^2+(a3+a4*u)*cos(u)+(a5+a6*u)*sin(u)]
                     % Out-of-plane: - G = G0 + G1*(u-u0) + G2*sin(u-u0+G3)
                     %               - H = H0 + H1*(u-u0) + H2*sin(u-u0+H3)
                     %               - h = h0 + h1*(u-u0) + h2*sin(u-u0+h3)
    
   if (param_flag == 0) % x = [a0 a1 a3 a5 G0 G1 H0 H1 h0 h1]
                        % l = [a2 a4 a6 G2 G3 H2 H3 h2 h3]
    
               parameters = [x(1) ; x(2) ; l(1) ; x(3) ; l(2) ; x(4) ; l(3); % a-parameters
                             x(5) ; x(6) ; l(4) ; l(5) ; % G-parameters
                             x(7) ; x(8) ; l(6) ; l(7) ; % H-parameters
                             x(9) ; x(10) ; l(8) ; l(9); % h-parameters
                             ];
                        
   end
   

elseif (shape_flag == 1) % Radial: 1/[a0+a1*u+a2*u^2+(a3+a4*u)*cos(u)+(a5+a6*u)*sin(u)]
                         % Out-of-plane: - G = G0 + G1*(u-u0) + G2*sin(u-u0+G3)
                         %               - i = i0 + i1*(u-u0) + i2*sin(u-u0+i3)
                         %               - h = h0 + h1*(u-u0) + h2*sin(u-u0+h3)
                         % Then, H = G*cos(i)
    
   if (param_flag == 0) % x = [a0 a1 a3 a5 G0 G1 i0 i1 h0 h1]
                        % l = [a2 a4 a6 G2 G3 i2 i3 h2 h3]
    
               parameters = [x(1) ; x(2) ; l(1) ; x(3) ; l(2) ; x(4) ; l(3); % a-parameters
                             x(5) ; x(6) ; l(4) ; l(5) ; % G-parameters
                             ];
                        
   elseif (param_flag == 1) % x = [a3 a4 a5 a6 G0 G1 i0 i1 h0 h1]
                            % l = [a0 a1 a2 G2 G3 i2 i3 h2 h3]
       
        parameters = [l(1) ; l(2) ; l(3) ; x(1) ; x(2) ; x(3) ; x(4); % a-parameters
                             x(5) ; x(6) ; l(4) ; l(5) ; % G-parameters
                             ];
   end
   
elseif (shape_flag == 2)
                 % a0   a1   a2   e0   e1   e2   w0   w1   w2   G0   G1   G2   G3 
    parameters = [x(1) x(2) l(1) x(3) x(4) l(2) x(5) x(6) l(3) x(7) x(8) l(4) l(5)]';
   
end

end