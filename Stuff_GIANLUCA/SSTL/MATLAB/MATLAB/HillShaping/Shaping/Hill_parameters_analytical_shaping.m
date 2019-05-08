function [r,G,v_r,H,h,u_dot,du_dot_du,i] = Hill_parameters_analytical_shaping(p,u0,shape_flag)



if (shape_flag == 0) % Radial: Inv-Polynomial ; Out-of-plane: Lin-Trig
    if nargout<8
    a0= p(1); a1= p(2); a2= p(3); a3 = p(4); a4 = p(5); a5 = p(6); a6=p(7);
    G0 = p(8); G1 = p(9); G2 = p(10); G3 = p(11);
    H0 = p(12); H1 = p(13); H2 = p(14); H3 = p(15);
    h0 = p(16); h1 = p(17); h2 = p(18); h3 = p(19);
   
    r = @(u) 1./( a0 + a1.*u + a2.*u.^2 + (a3+a4.*u).*cos(u) + (a5+a6.*u).*sin(u) );
    
    G = @(u) G0 + G1.*(u-u0) + G2.*sin(u-u0+G3);
    
    H = @(u) H0 + H1.*(u-u0) + H2.*sin(u-u0+H3);
    
    h = @(u) h0 + h1.*(u-u0) + h2.*cos(u-u0+h3);
    
    v_r = @(u) -0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u));
    
    syms f_z
    % Real
%     u_dot = @(u) (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 - 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .* sin(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.1e1 ./ 0.2e1);
%     du_dot_du = @(u) (G1 + G2 .* cos(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 + 0.2e1 .* (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .* sin(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.1e1 ./ 0.2e1) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) - 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (H1 + H2 .* cos(u - u0 + H3)) .* sin(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.1e1 ./ 0.2e1) - 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .* cos(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.1e1 ./ 0.2e1) + 0.2e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .* sin(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 3 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.1e1 ./ 0.2e1) .* (G1 + G2 .* cos(u - u0 + G3)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .* sin(u) .* f_z ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (0.1e1 - (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2) .^ (-0.3e1 ./ 0.2e1) .* (-0.2e1 .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 .* (H1 + H2 .* cos(u - u0 + H3)) + 0.2e1 .* (H0 + H1 .* (u - u0) + H2 .* sin(u - u0 + H3)) .^ 2 ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 3 .* (G1 + G2 .* cos(u - u0 + G3))) ./ 0.2e1;

    % Approximated
    u_dot = @(u) (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2;
    du_dot_du = @(u) (G1 + G2 .* cos(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 + 0.2e1 .* (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u));
    end
    
elseif ( shape_flag == 1 )
    
    a0= p(1); a1= p(2); a2= p(3); a3 = p(4); a4 = p(5); a5 = p(6); a6=p(7);
    G0 = p(8); G1 = p(9); G2 = p(10); G3 = p(11);
    if length(p) == 19
        i0 = p(12); i1 = p(13); i2 = p(14); i3 = p(15);
        h0 = p(16); h1 = p(17); h2 = p(18); h3 = p(19);
    end
        
    
    r = @(u) 1./( a0 + a1.*u + a2.*u.^2 + (a3+a4.*u).*cos(u) + (a5+a6.*u).*sin(u) );
    
    G = @(u) G0 + G1.*(u-u0) + G2.*sin(u-u0+G3);
    
    if nargout > 2
    i = @(u) i0 + i1.*(u-u0) + i2.*sin(u-u0+i3);
    
    H = @(u) (G0 + G1.*(u-u0) + G2.*sin(u-u0+G3)).*cos(i0 + i1.*(u-u0) + i2.*sin(u-u0+i3));
    
    h = @(u) h0 + h1.*(u-u0) + h2.*cos(u-u0+h3);
    
    v_r = @(u) -0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u));
    
    syms f_z
    % Real % Qui va cambiato con nuova parametrizzazione
%     u_dot = @(u) (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 - 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* cos(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* sin(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) ./ sin(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3));
%     du_dot_du = @(u) (G1 + G2 .* cos(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 + 0.2e1 .* (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 .* cos(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* sin(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) ./ sin(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (i1 + i2 .* cos(u - u0 + i3)) .* sin(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) - 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* cos(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* cos(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) ./ sin(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* cos(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* sin(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .^ 2 ./ sin(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .* (G1 + G2 .* cos(u - u0 + G3)) + 0.1e1 ./ (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* cos(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .^ 2 .* sin(u) .* f_n ./ (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) ./ sin(i0 + i1 .* (u - u0) + i2 .* sin(u - u0 + i3)) .^ 2 .* (i1 + i2 .* cos(u - u0 + i3));

    % Approximated
    u_dot = @(u) (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 ;

    du_dot_du = @(u) (G1 + G2 .* cos(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .^ 2 + 0.2e1 .* (G0 + G1 .* (u - u0) + G2 .* sin(u - u0 + G3)) .* (a0 + a1 .* u + a2 .* u .^ 2 + (a4 .* u + a3) .* cos(u) + (a6 .* u + a5) .* sin(u)) .* (a1 + 0.2e1 .* a2 .* u + a4 .* cos(u) - (a4 .* u + a3) .* sin(u) + a6 .* sin(u) + (a6 .* u + a5) .* cos(u));

    end
end


end