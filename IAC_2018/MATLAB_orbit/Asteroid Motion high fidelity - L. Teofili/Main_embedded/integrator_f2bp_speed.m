function dX_dt = integrator_f2bp_speed_copy(t,p,varargin)

global r1b0 r2b0 r3b0 r1a0 r2a0 r3a0 Ja Jb rho_a rho_b nfaces_a nfaces_b I1 I2 norm_mass_a norm_mass_b c1a c2a c3a c1b c2b c3b;

% ______ renaming the variables

rcm_a = p(1:3);
vcm_a = p(4:6);
rcm_b = p(7:9);
vcm_b = p(10:12);
w_a = p(13:15);
w_b = p(16:18);
En_tot = p(37);


S = [p(19:21)';p(22:24)';p(25:27)'];
P = [p(28:30)';p(31:33)';p(34:36)'];

%______________Computation

a_p = P(:,1);
b_p = P(:,2);
c_p = P(:,3);

St = S';
Pt = P';

r1b = S*r1b0;
r2b = S*r2b0;
r3b = S*r3b0;

r1a = P*r1a0;
r2a = P*r2a0;
r3a = P*r3a0;

R = (rcm_b-rcm_a);             % Vector starting on A ending in B

%% Testing Function
%[Mb_T,Mb_f_T,E_T] = torque_action_calculator(R,S,P,r1a0,r2a0,r3a0,r1b0,r2b0,r3b0,Ja,Jb,rho_a,rho_b,nfaces_a,nfaces_b)
%[Fanp,Enp,Unp] = Action_calculator_forces2_torques2_ok_pluspotential(R,r1a,r2a,r3a,r1a0,r2a0,r3a0,r1b,r2b,r3b,Ja,Jb,rho_a,rho_b,nfaces_a,nfaces_b);

%% Tested Function
tic
% [Fa,E,U] = Action_calculator_forces2_torques2_ok_pluspotential_par(R,r1a,r2a,r3a,r1a0,r2a0,r3a0,r1b,r2b,r3b,Ja,Jb,rho_a,rho_b,nfaces_a,nfaces_b);
[Fa,E,U] = Action_calculator_forces2_torques2_ok_pluspotential(R,r1a,r2a,r3a,r1a0,r2a0,r3a0,r1b,r2b,r3b,Ja,Jb,rho_a,rho_b,nfaces_a,nfaces_b);
toc

%% Force and Torques

Fb = -Fa;
Ma = (cross(a_p,E(:,1))+cross(b_p,E(:,2))+cross(c_p,E(:,3)));     %  torque acting on the A body expressed in the inertial frame refer to my
                                                                  %  article not the Fanhestock one 
Mb = - Ma - cross(R,Fb);
Ma = Pt*Ma;
Mb = St*Mb;

%% Forces and torques check

val_conf_f = max(abs(Fb));
Fb(abs(Fb)./val_conf_f <=1e-5)=0;
% Fb(abs(Fb)<1e-10)=0;
val_conf_f = max(abs(Fa));
Fa(abs(Fa)./val_conf_f <=1e-5)=0;
% Fa(abs(Fa)<1e-10)=0;
val_conf = max(abs(Mb));
Mb(abs(Mb)./val_conf <=1e-6)=0;
% Mb(abs(Ma)<1e-10)=0;
val_conf = max(abs(Ma));
Ma(abs(Ma)./val_conf <=1e-6)=0;
% Ma(abs(Ma)<1e-10)=0;


%% Energy box
norms_vcm_a = vcm_a(1)^2+vcm_a(2)^2+vcm_a(3)^2;
norms_vcm_b = vcm_b(1)^2+vcm_b(2)^2+vcm_b(3)^2;
w_a_in = P*w_a;
w_b_in = S*w_b;
I1_in = P*I1*Pt;
I2_in = S*I2*St;

En_tot = -U+0.5*norm_mass_a*norms_vcm_a+0.5*norm_mass_b*norms_vcm_b+0.5*w_a_in'*I1_in*w_a_in+0.5*w_b_in'*I2_in*w_b_in;
En_cin = 0.5*w_a_in'*I1_in*w_a_in+0.5*w_b_in'*I2_in*w_b_in+0.5*norm_mass_a*norms_vcm_a+0.5*norm_mass_b*norms_vcm_b;


%% Computing acceleration

Mb_I1=Mb(1)/I2(1,1);
Mb_I2=Mb(2)/I2(2,2);
Mb_I3=Mb(3)/I2(3,3);

Ma_I1=Ma(1)/I1(1,1);
Ma_I2=Ma(2)/I1(2,2);
Ma_I3=Ma(3)/I1(3,3);


w_p_b_x = c1b*w_b(2)*w_b(3)+ Mb_I1;
w_p_b_y = c2b*w_b(1)*w_b(3)+ Mb_I2;
w_p_b_z = c3b*w_b(2)*w_b(1)+ Mb_I3;

w_p_b = [w_p_b_x;w_p_b_y;w_p_b_z];
w_p_b(abs(w_p_b) <=1e-5)=0;

w_p_a_x = c1a*w_a(2)*w_a(3)+ Ma_I1;
w_p_a_y = c2a*w_a(1)*w_a(3)+ Ma_I2;
w_p_a_z = c3a*w_a(2)*w_a(1)+ Ma_I3;

w_p_a = [w_p_a_x;w_p_a_y;w_p_a_z];
w_p_a(abs(w_p_a) <=1e-5)=0;


wa_x_=skw(w_a);
wb_x_=skw(w_b);

dS_dT = S*wb_x_;
dP_dT = P*wa_x_;

t/(2*pi)

dX_dt=[
    vcm_a;
    Fa/norm_mass_a;
    vcm_b;
    Fb/norm_mass_b;
    w_p_a;
    w_p_b;
    dS_dT(1,1);
    dS_dT(1,2);
    dS_dT(1,3);
    dS_dT(2,1);
    dS_dT(2,2);
    dS_dT(2,3);
    dS_dT(3,1);
    dS_dT(3,2);
    dS_dT(3,3);
    dP_dT(1,1);
    dP_dT(1,2);
    dP_dT(1,3);
    dP_dT(2,1);
    dP_dT(2,2);
    dP_dT(2,3);
    dP_dT(3,1);
    dP_dT(3,2);
    dP_dT(3,3);
    En_tot;
    -U;
    En_cin;
    ];

end


















