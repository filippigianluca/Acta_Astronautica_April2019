%Main_Asteroid_Motion Computes the motion of two asteroids of whatever shape
%  one around the other given the initial conditions. The code is based on
%  the algorithm exposed in "Periodic Orbits in the Binary Asteroid System 1999KW4"
%  Lorenzo Teofili et al., a further reference is "Simulation of the full
%  two rigid body problem using polyhedral mutual potential and potential
%  derivatives approach" E. G. Fahnestock, D. J. Scheeres.
%  Lorenzo Teofili, University of Rome La Sapienza, University of Strathclyde, 5/11/2015

clear all
close all
clc

global G r1a r2a r3a r1b0 r2b0 r3b0 r1a0 r2a0 r3a0 Ja Jb rho_a rho_b nfaces_a nfaces_b Pt S a_s b_s c_s w_p_b_2 w_p_a_1 i t_wp I1 I2 norm_mass_a norm_mass_b c1a c2a c3a c1b c2b c3b ni;

%% Load Asteroid Model
% Load asteroid model from database data
% use ellipsoid_mesh(a_x,a_y,a_z) to load an ellipses of semimajor axes a_x,a_y and a_z
% use sphere_mesh(r) to load a sphere with radious r 

[r11,r21,r31,cent_vert_a,faces_a] = load_model_2bp('data/kw4a.tab');%load_cube_2bp(14,1,1);%ellipsoid_mesh(5,0.1,0.1);%load_model_2bp('data/kw4a.tab');%ellipsoid_mesh(2,1,1);%sphere_mesh(0.04);%load_model_2bp('data/kw4a.tab');%%%load_model_cylinder([],[]);%ellipsoid_mesh(14,0.1,0.1)%sphere_mesh(1);%load_model_2bp('data/kw4b.tab');%ellipsoid_mesh(0.285,0.2225,0.1715);%load_model_2bp('data/kw4b.tab');%load_model_2bp('data/kw4b.tab');%ellipsoid_mesh(0.4,0.1,0.1);%ellipsoid_mesh(0.685,0.2225,0.1715);%ellipsoid_mesh(0.685,0.2225,0.1715);%sphere_mesh(0.768);%ellipsoid_mesh(0.285,0.2225,0.1715);%load_model_2bp('data/kw4b.tab');%sphere_mesh(1738,R_cm2)%
[r12,r22,r32,cent_vert_b,faces_b] = load_model_2bp('data/kw4b.tab');%load_cube_2bp(0.001,0.001,0.001);%load_model_2bp('data/kw4b.tab');%ellipsoid_mesh(5,1,1);%load_model_2bp('data/kw4b.tab');%sphere_mesh(1);%load_model_2bp('data/kw4a.tab');%sphere_mesh(0.1);%load_model_2bp('data/kw4a.tab');%sphere_mesh(0.468);%sphere_mesh(0.768);%load_model_2bp('data/kw4a.tab',R_cm2)%sphere_mesh(6378,R_cm1)%

%% ---------------Defining Independent Variables---------------%%

G = 6.67384 * 1e-20;             % Gravitational constant
deg2rad = pi/180;                % Convertion factor from deg to radiant

mass_a = 2.353*10^(12);          % kg mass of the body A
mass_b = 0.135*10^(12);          % kg mass of the body B
distance = 2.548;                % Km distance between the asteroid centres of mass
direction=[1,0,0];               % Relative position unit vector of the asteroid centres of mass

theta_in1 = deg2rad*[0,0,0];     % Initial attitude of body A in euler angles [x_angle,y_angle,z_angle] 
theta_in2 = deg2rad*[0,-180,0];  % Initial attitude of body B in euler angles [x_angle,y_angle,z_angle]
w_in1_rot = [0,0,0];             % rad/s Initial angular velocity in the inertial reference frame of body A
w_in2_rot = [0,0,0];             % rad/s Initial angular velocity in the inertial reference frame of body B

circular_motion = 1;             % Put equal to 1 to set automatically the condition of circular motion
V_cm1 = [0, 0, 0];               % Km/s Barycentre velocity of the body A
V_cm2 = [0, 0, 0];               % Km/s Barycentre velocity of the body B

t_f = 22*2*pi;                   % Final Time for the integration 2*pi correnspond to a period 
                                 % of the primary around the secondary bodies in an ipotetical round orbit

RelTol = 1e-13;                  % RelTol in ODE45
AbsTol = 1e-13;                  % AbsTol in ODE45

%% Defining dependent Variables

%___ adimensional Variables

norm_mass = (mass_a+mass_b);
norm_dist = distance;
norm_time = sqrt(distance^3/(G*norm_mass));
norm_velocity = norm_dist/norm_time;
norm_ang_velocity = 1/norm_time;
norm_density = norm_mass/norm_dist^3;

r11 = r11/norm_dist;
r21 = r21/norm_dist;
r31 = r31/norm_dist;
r12 = r12/norm_dist;
r22 = r22/norm_dist;
r32 = r32/norm_dist;

norm_mass_a = mass_a/norm_mass;
norm_mass_b = mass_b/norm_mass;

R = direction/norm(direction);

V_cm1 = V_cm1/norm_velocity;
V_cm2 = V_cm2/norm_velocity;
w_in1_rot = w_in1_rot/norm_ang_velocity;
w_in2_rot = w_in2_rot/norm_ang_velocity;

sum_mass = mass_a+mass_b;

ni = mass_a/(mass_a+mass_b);

R_cm_a = -(1-ni)*R;   % Position of the center of mass 1 respect to an inertial reference frame (Usefull just for representation not utilized )
R_cm_b = (ni)*R;      % Position of the center of mass 2 respect to an inertial reference frame (Usefull just for representation)

[I,volume_a,Ja]=Asteroid_Inertia_matrix_2(r11,r21,r31);
[I,volume_b,Jb]=Asteroid_Inertia_matrix_2(r12,r22,r32);

rho_a = norm_mass_a/volume_a;
rho_b = norm_mass_b/volume_b;

nfaces_a=length(Ja);
nfaces_b=length(Jb);

[I1,volume,mass1]=Asteroid_Inertia_matrix(r11,r21,r31,rho_a);
[I2,volume,mass2]=Asteroid_Inertia_matrix(r12,r22,r32,rho_b);

I1 = diag(diag(I1));
I2 = diag(diag(I2));

c1a = -(I1(3,3)-I1(2,2))/I1(1,1);
c2a = -(I1(1,1)-I1(3,3))/I1(2,2);
c3a = -(I1(2,2)-I1(1,1))/I1(3,3);
c1b = -(I2(3,3)-I2(2,2))/I2(1,1);
c2b = -(I2(1,1)-I2(3,3))/I2(2,2);
c3b = -(I2(2,2)-I2(1,1))/I2(3,3);

r1b0 =r12';
r2b0 =r22';
r3b0 =r32';

r1a0 =r11';
r2a0 =r21';
r3a0 =r31';


%% Initial conditions for the Integration

R_cm1 = R_cm_a;
R_cm2 = R_cm_b;

if circular_motion == 1;
    V_cm1_circ_rel = sqrt(G*sum_mass/distance);
    w = V_cm1_circ_rel*norm_time/distance;
    V_circ_moon = w*norm(R_cm_b);
    V_circ_earth = w*norm(R_cm_a);
    V_cm1 = [0, -V_circ_earth, 0];
    V_cm2 = [0, V_circ_moon, 0];
    theta_in1 = [0*pi/180,0*pi/180,0*pi/180]; 
    theta_in2 = [0*pi/180,-180*pi/180,0*pi/180];
    walpha = 6.301543631863305e+00;
    wbeta = norm(V_cm1)/norm(R_cm1);
    w_in1_rot = [0,0,walpha];             
    w_in2_rot = [0,0,-wbeta];             
end

P = rotz(theta_in1(3)*180/pi)*roty(theta_in1(2)*180/pi)*rotx(theta_in1(1)*180/pi);
S = rotz(theta_in2(3)*180/pi)*roty(theta_in2(2)*180/pi)*rotx(theta_in2(1)*180/pi);
Pt = P';
St = S';

vers_V_cm1 = 0.5*V_cm1/norm(V_cm1);
vers_V_cm2 = 0.5*V_cm2/norm(V_cm2);


%%  Plot The Situation
%__ Plot the Initial situation in the reference frame of the A body

% Initial_attitude_a = (P'*S*[cent_vert_a(:,1), cent_vert_a(:,2),cent_vert_a(:,3)]')';
% Initial_attitude_b = ([cent_vert_b(:,1), cent_vert_b(:,2),cent_vert_b(:,3)]')';
% initial_R_cm_a     =  P'*S*R_cm_a';
% initial_R_cm_b     =  P'*S*R_cm_b';

%__ Plot the Initial situation in the inertial reference frame

Initial_attitude_a = ((P*[cent_vert_a(:,1), cent_vert_a(:,2),cent_vert_a(:,3)]')')./norm_dist;
Initial_attitude_b = ((S*[cent_vert_b(:,1), cent_vert_b(:,2),cent_vert_b(:,3)]')')./norm_dist;
initial_R_cm_a     =  R_cm_a';
initial_R_cm_b     =  R_cm_b';

fva.vertices = [Initial_attitude_a(:,1) + initial_R_cm_a(1), Initial_attitude_a(:,2) + initial_R_cm_a(2),Initial_attitude_a(:,3) + initial_R_cm_a(3)];
fva.faces = faces_a;
fva.facevertexcdata = zeros(nfaces_a,1);
fvb.vertices = [Initial_attitude_b(:,1) + initial_R_cm_b(1), Initial_attitude_b(:,2) + initial_R_cm_b(2),Initial_attitude_b(:,3) + initial_R_cm_b(3)];
fvb.faces = faces_b;
fvb.facevertexcdata = zeros(nfaces_b,1);


%plot both asteroids
figure(300)
title('Inertial Reference Frame')
patch(fva,'FaceColor','flat')
patch(fvb,'FaceColor','flat')
axis equal
axis([-inf,inf,-inf,inf,-inf,inf])

% Plot asteroids name
% text(max(0*Initial_attitude_a(:,1) + initial_R_cm_a(1)),max(0*Initial_attitude_a(:,2) + initial_R_cm_a(2)),max(1.2*Initial_attitude_a(:,3) + initial_R_cm_a(3)),'\bf{KW\alpha}');
% text(max(0*Initial_attitude_b(:,1) + initial_R_cm_b(1)),max(0*Initial_attitude_b(:,2) + initial_R_cm_b(2)),max(1.4*Initial_attitude_b(:,3) + initial_R_cm_b(3)),'\bf{KW\beta}');

%%Plot velocities
hold on
% quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),vers_V_cm1(1),vers_V_cm1(2),vers_V_cm1(3),'k'); text(initial_R_cm_a(1)+vers_V_cm1(1),initial_R_cm_a(2)+vers_V_cm1(2),initial_R_cm_a(3)+vers_V_cm1(3),'\bf{V_{KW\alpha}}')
% quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),vers_V_cm2(1),vers_V_cm2(2),vers_V_cm2(3),'k'); text(initial_R_cm_b(1)+vers_V_cm2(1),initial_R_cm_b(2)+vers_V_cm2(2),initial_R_cm_b(3)+vers_V_cm2(3),'\bf{V_{KW\beta}}')

xlabel('x')
ylabel('y')
zlabel('z')
%grid

%Plot Attitude A
e1 = [1,0,0];
e2 = [0,1,0];
e3 = [0,0,1];
e1r = 0.7*(P*e1')';
e2r = 0.7*(P*e2')';
e3r = 0.7*(P*e3')';
figure(301)
title('A Body Attitude')
patch(fva,'FaceColor','white')
hold on
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e1(1),e1(2),e1(3),'k'); text(initial_R_cm_a(1)+e1(1),initial_R_cm_a(2)+e1(2),initial_R_cm_a(3)+e1(3),'X')
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e2(1),e2(2),e2(3),'k'); text(initial_R_cm_a(1)+e2(1),initial_R_cm_a(2)+e2(2),initial_R_cm_a(3)+e2(3),'Y')
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e3(1),e3(2),e3(3),'k'); text(initial_R_cm_a(1)+e3(1),initial_R_cm_a(2)+e3(2),initial_R_cm_a(3)+e3(3),'Z')
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e1r(1),e1r(2),e1r(3),'--k'); text(initial_R_cm_a(1)+e1r(1),initial_R_cm_a(2)+e1r(2),initial_R_cm_a(3)+e1r(3),'x')
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e2r(1),e2r(2),e2r(3),'--k'); text(initial_R_cm_a(1)+e2r(1),initial_R_cm_a(2)+e2r(2),initial_R_cm_a(3)+e2r(3),'y')
quiver3(initial_R_cm_a(1),initial_R_cm_a(2),initial_R_cm_a(3),e3r(1),e3r(2),e3r(3),'--k'); text(initial_R_cm_a(1)+e3r(1),initial_R_cm_a(2)+e3r(2),initial_R_cm_a(3)+e3r(3),'z')
axis equal

%Plot Attitude B
e1 = 0.5*[1,0,0];
e2 = 0.5*[0,1,0];
e3 = 0.5*[0,0,1];
e1r = 0.4*(S*e1')';
e2r = 0.4*(S*e2')';
e3r = 0.4*(S*e3')';
figure(302)
title('B Body Attitude')
patch(fvb,'FaceColor','white')
hold on
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e1(1),e1(2),e1(3),'k'); text(initial_R_cm_b(1)+e1(1),initial_R_cm_b(2)+e1(2),initial_R_cm_b(3)+e1(3),'X')
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e2(1),e2(2),e2(3),'k'); text(initial_R_cm_b(1)+e2(1),initial_R_cm_b(2)+e2(2),initial_R_cm_b(3)+e2(3),'Y')
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e3(1),e3(2),e3(3),'k'); text(initial_R_cm_b(1)+e3(1),initial_R_cm_b(2)+e3(2),initial_R_cm_b(3)+e3(3),'Z')
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e1r(1),e1r(2),e1r(3),'--k'); text(initial_R_cm_b(1)+e1r(1),initial_R_cm_b(2)+e1r(2),initial_R_cm_b(3)+e1r(3),'x')
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e2r(1),e2r(2),e2r(3),'--k'); text(initial_R_cm_b(1)+e2r(1),initial_R_cm_b(2)+e2r(2),initial_R_cm_b(3)+e2r(3),'y')
quiver3(initial_R_cm_b(1),initial_R_cm_b(2),initial_R_cm_b(3),e3r(1),e3r(2),e3r(3),'--k'); text(initial_R_cm_b(1)+e3r(1),initial_R_cm_b(2)+e3r(2),initial_R_cm_b(3)+e3r(3),'z')
axis equal



%% Integration of motion

tspan = [0, t_f];
in_condition = [R_cm1,V_cm1,R_cm2,V_cm2,w_in1_rot,w_in2_rot,S(1,:),S(2,:),S(3,:),P(1,:),P(2,:),P(3,:),0,0,0];
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
%[tsp,p]=ode45('integrator_f2bp_ok2',tspan,in_condition,options)%,I1,I2,mass_a,mass_b,c1a,c2a,c3a,c1b,c2b,c3b);
sol = ode113(@integrator_f2bp_speed,tspan,in_condition,options);
t = 0:0.01:t_f;
[x,xp] = deval(sol,t);

%% Save data

%  save('asteroid_variables','p','cent_vert_a','cent_vert_b','faces_a','faces_b','tsp','norm_dist');
%  ast_movie_in_cond

%% Postprocess Data


tsp = 0:1:(length(x(1,:))-1);
tsp = tsp*t_f/(length(x(1,:))-1);
norm_Ra = sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2);
norm_Rb = sqrt(x(7,:).^2+x(8,:).^2+x(9,:).^2);
norm_Va = sqrt(xp(1,:).^2+xp(2,:).^2+xp(3,:).^2);
norm_Vb = sqrt(xp(7,:).^2+xp(8,:).^2+xp(9,:).^2);
pos_com = (mass_a*[x(1,:);x(2,:);x(3,:)]+mass_b*[x(7,:);x(8,:);x(9,:)])./(mass_a+mass_b);
C11P = x(28,:);
C21P = x(29,:);
C31P = x(30,:);
C32P = x(33,:);
C33P = x(36,:);
[phiA,theA,psiA]= att2eulang(C11P,C21P,C31P,C32P,C33P);
C11S = x(19,:);
C21S = x(20,:);
C31S = x(21,:);
C32S = x(24,:);
C33S = x(27,:);
[phiB,theB,psiB]= att2eulang(C11S,C21S,C31S,C32S,C33S);
phiB(phiB<0)=phiB(phiB<0)+2*pi;


% %% PLOT
%
figure(2)
plot3(x(1,:),x(2,:),x(3,:))
hold on
plot3(x(7,:),x(8,:),x(9,:),'r')
plot3(pos_com(:,1),pos_com(:,2),pos_com(:,3),'+')
% plot3(pos_com(:,1),pos_com(:,2),pos_com(:,3))
% plot3(R_cm1(1),R_cm1(2),R_cm1(3),'o')
% plot3(R_cm2(1),R_cm2(2),R_cm2(3),'+')
xlabel('X')
ylabel('Y')
%zlabel('Z')
axis equal
legend('Body A Trajectory','Body B Trajectory','Center of Mass')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1traject
%
figure(3)
plot(tsp,x(13,:))
hold on
plot(tsp,x(14,:))
plot(tsp,x(15,:))
legend('\omega_{Ax}','\omega_{Ay}','\omega_{Az}')
%title('A Angular Velocity on Body Axes')
xlabel('Time [non dimensional]')
ylabel('Angular Velocity Components [non dimensional]')
grid
% print -depsc /Users/Lorenzo/Desktop/figure_fragart/1wa
%
figure(4)
plot(tsp,x(16,:))
hold on
plot(tsp,x(17,:))
plot(tsp,x(18,:))
legend('\omega_{Bx}','\omega_{By}','\omega_{Bz}')
%title('B Angular Velocity on Body Axes')
xlabel('Time [non dimensional]')
ylabel('Angular Velocity Components [non dimensional]')
grid
%print -depsc /Users/Lorenzo/Desktop/figure_fragart/1wb
%
figure(5)
plot(tsp,180/pi*psiA)
hold on
plot(tsp,180/pi*theA)
plot(tsp,180/pi*phiA)
legend('\psi_{A}','\theta_{A}','\phi_{A}')
%title('A Euler Angles')
xlabel('Time [non dimensional]')
ylabel('Euler Angles [deg]')
grid
% print -depsc /Users/Lorenzo/Desktop/figure_fragart/1eulA
%
figure(6)
plot(tsp,180/pi*psiB)
hold on
plot(tsp,180/pi*theB)
plot(tsp,180/pi*phiB)
legend('\psi_{B}','\theta_{B}','\phi_{B}')
%title('B Euler Angles')
xlabel('Time [non dimensional]')
ylabel('Euler Angles [deg]')
grid
% print -depsc /Users/Lorenzo/Desktop/figure_fragart/1eulb
%
figure(7)
plot(tsp,15*norm_Ra)
hold on
plot(tsp,norm_Rb)
legend('r_A','r_B')
xlabel('Time [non dimensional]')
ylabel('\bf{r} \rm{Modulus [non dimensional]}')
%title('Radious of the Orbit')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1normR
%
figure(8)
plot(tsp,xp(1,:))
hold on
plot(tsp,xp(2,:))
plot(tsp,xp(3,:))
legend('v_{Ax}','v_{Ay}','v_{Az}')
%title('A Inertial Velocities Velocities ')
xlabel('Time [non dimensional]')
ylabel('Velocity [non dimensional]')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1compVa
%
figure(9)
plot(tsp,xp(7,:))
hold on
plot(tsp,xp(8,:))
plot(tsp,xp(9,:))
legend('v_{Bx}','v_{By}','v_{Bz}')
%title('B Inertial Velocities Velocities ')
xlabel('Time [non dimensional]')
ylabel('Velocity [non dimensional]')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1compVb
%
figure(10)
plot(tsp,15*norm_Va)
hold on
plot(tsp,norm_Vb)
legend('v_A','v_B')
%title('Velocity along the Orbit')
xlabel('Time [non dimensional]')
ylabel('\bf{v} \rm{Modulus [non dimensional]}')
%ylabel(['\bf{first line} \rm{second line}'])
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1normV

figure(11)
plot(tsp,xp(end-2,:))
hold on
plot(tsp,xp(end-1,:))
plot(tsp,xp(end,:))
legend('Total Energy','Potential Energy','Kinetic Energy')
%title('Energies of the System ')
xlabel('Time [non dimentional]')
ylabel('Energy [non dimensional]')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1energy
%
figure(12)
plot(tsp,x(7,:))
hold on
plot(tsp,x(8,:))
plot(tsp,x(9,:))
legend('X_{B}','Y_{B}','Z_{B}')
%title('B Inertial Position ')
xlabel('Time [non dimensional]')
ylabel('Position [non dimensional]')
grid
%  print -depsc /Users/Lorenzo/Desktop/figure_fragart/1posB
%
figure(13)
plot(tsp,x(1,:))
hold on
plot(tsp,x(2,:))
plot(tsp,x(3,:))
legend('X_{A}','Y_{A}','Z_{A}')
%title('A Inertial Position  ')
xlabel('Time [non dimensional]')
ylabel('Position [non dimensional]')
grid
