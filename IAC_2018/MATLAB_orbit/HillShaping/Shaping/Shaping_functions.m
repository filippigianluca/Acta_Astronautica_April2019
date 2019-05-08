function [c,A,b] = Shaping_functions(dep,arr,shape_flag,param_flag,lam)

% Shaping_functions: Function that construct the A-matrix, vector b and c
%                    that describes the boundary relations as:
%                    A*x + c(l) = b
%                    Where x is the vector with the shaping parameters to
%                    be found.
%
% INPUTs:  - Hill_dep: Hills parameters at departure [6x1]
%          - Hill_arr: Hills parameters at arrival [6x1]
%          - shape_flag: Flag to indicate which shape is selected [-]
%          - param_flag: Flag indicating which parameters are used for BC
%                        and which for optimization
%          - l: Lambda, vector with values of free parameters for optimization
%          - nr: Number of revolutions
%
% OUTPUTs: - c: Vector describing the relation in the boundary relations
%               which depend only on l parameters (to be optimized and 
%               in this function supposed given)
%          - A: Matrix describing the linear relation of the BCs
%               with the shaping parameters
%          - b: Known vector with the Hill parameters at BCs
%
% N.B. c-vector is the first output as when l changes A & b remain the
%      same and we do not need to compute them again (nargout will be used).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by: Cristian Greco, 08-09-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Departure argument of latitude
u_0 = dep(4);
cos_u0 = cos(u_0);
sin_u0 = sin(u_0);

% Arrival argument of latitude
u_f = arr(4);
cos_uf = cos(u_f);
sin_uf = sin(u_f);

du = u_f - u_0;

if shape_flag == 1
    
    if param_flag == 0
    
    c = [ lam(1).*u_0.^2 + lam(2).*u_0.*cos_u0 + lam(3).*u_0.*sin_u0 ; % 1/r0
        lam(1).*u_f.^2 + lam(2).*u_f.*cos_uf + lam(3).*u_f.*sin_uf ; % 1/rf
        2*lam(1)*u_0 + lam(2).*(-u_0.*sin_u0+cos_u0) + lam(3).*(u_0.*cos_u0+sin_u0); % -v_r0/r0^2*u0_dot
        2*lam(1)*u_f + lam(2).*(-u_f.*sin_uf+cos_uf) + lam(3).*(u_f.*cos_uf+sin_uf); % -v_rf/rf^2*uf_dot
        lam(4).*sin(lam(5)) ; % G(t0)
        lam(4).*sin(du+lam(5)); % G(tf)
        ];
    
    A = [1 u_0 cos_u0 sin_u0 0 0 ;
        1 u_f cos_uf sin_uf 0 0 ;
        0 1 -sin_u0 cos_u0 0 0 ;
        0 1 -sin_uf cos_uf 0 0 ;
        0 0 0 0 1 0 ;
        0 0 0 0 1 du];
    
    u0_dot = dep(3)/dep(1)^2; %theta_dot = G/(r^2)=Angular momentum/r^2 (neglecting omega_dot)
    uf_dot = arr(3)/arr(1)^2;
    
    b = [1/dep(1); % 1/r0
        1/arr(1); % 1/rf
        -dep(2)/dep(1)^2/u0_dot; % -v_r0/r0^2*u0_dot
        -arr(2)/arr(1)^2/uf_dot; % -v_rf/rf^2*uf_dot
        dep(3); % G(t0)
        arr(3); % G(tf)
        ];
    

    elseif param_flag == 1
        
        c = [ lam(1) + lam(2).*u_0 + lam(3).*u_0.^2 ; % 1/r0
        lam(1) + lam(2).*u_f + lam(3).*u_f.^2 ; % 1/rf
        lam(2) + 2*lam(3).*(u_0); % -v_r0/r0^2*u0_dot
        lam(2) + 2*lam(3).*(u_f); % -v_rf/rf^2*uf_dot        
        lam(4).*sin(lam(5)) ; % G(t0)
        lam(4).*sin(du+lam(5)); % G(tf)
        ];
    
    
    A = [cos_u0 u_0*cos_u0 sin_u0 u_0*sin_u0 0 0 ;
        cos_uf u_f*cos_uf sin_uf u_f*sin_uf 0 0 ;
        -sin_u0 (cos_u0-u_0*sin_u0) cos_u0 (sin_u0+u_0*cos_u0) 0 0 ;
        -sin_uf (cos_uf-u_f*sin_uf) cos_uf (sin_uf+u_f*cos_u0) 0 0 ;
        0 0 0 0 1 0 ;
        0 0 0 0 1 du];
    
    u0_dot = dep(3)/dep(1)^2; %theta_dot = G/(r^2)=Angular momentum/r^2 (neglecting omega_dot)
    uf_dot = arr(3)/arr(1)^2;
    
    b = [1/dep(1); % 1/r0
        1/arr(1); % 1/rf
        -dep(2)/dep(1)^2/u0_dot; % -v_r0/r0^2*u0_dot
        -arr(2)/arr(1)^2/uf_dot; % -v_rf/rf^2*uf_dot
        dep(3); % G(t0)
        arr(3); % G(tf)
        ];
%     keyboard
        
    end
    
elseif shape_flag == 2
    
    
    c = [ 0 ;
          lam(1)*sin(du);
          0;
          lam(2)*sin(du);
          0;
          lam(3)*sin(du);
          lam(4).*sin(lam(5)) ; % G(t0)
          lam(4).*sin(du+lam(5)); % G(tf)
        ];
    
    A = [1 0 0 0 0 0 0 0 ; % a(t0)
        1 du 0 0 0 0 0 0; % a(tf)
        0 0 1 0 0 0 0 0 ; % e(t0)
        0 0 1 du 0 0 0 0; % e(tf)
        0 0 0 0 1 0 0 0; % w(t0)
        0 0 0 0 1 du 0 0; % w(tf)
        0 0 0 0 0 0 1 0 ; % G(t0)
        0 0 0 0 0 0 1 du]; % G(tf)
    
    
    b = [dep(1) ; % a0
        arr(1); % af
        dep(2); % e0
        arr(2); % ef
        dep(3); % w0
        arr(3); % wf
        dep(4); % G0
        arr(4); % Gf
        ];
    
end

end

