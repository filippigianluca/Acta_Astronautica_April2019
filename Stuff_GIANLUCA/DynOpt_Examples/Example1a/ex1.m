clear all; close all; clc

s = 'http://byu.apmonitor.com';
a = 'ex1a_dynopt';

addpath('apm')
apm(s,a,'clear all');
apm_load(s,a,'ex1.apm');
csv_load(s,a,'ex1.csv');

apm_option(s,a,'nlc.nodes',4);
apm_option(s,a,'nlc.solver',1);
apm_option(s,a,'nlc.imode',6);
apm_option(s,a,'nlc.mv_type',1);

apm_info(s,a,'MV','u');
apm_option(s,a,'u.status',1);
apm_option(s,a,'u.dcost',0);

output = apm(s,a,'solve');
disp(output)
y = apm_sol(s,a); z = y.x;

disp(['Optimal Solution: ' num2str(z.x2(end))])

figure(1)

subplot(2,1,1)
plot(z.time,z.u,'r-','LineWidth',2)
legend('u')
ylabel('Manipulated')

subplot(2,1,2)
plot(z.time,z.x1,'b--','LineWidth',2)
hold on
plot(z.time,z.x2,'g:','LineWidth',2)
legend('x_1','x_2')
ylabel('Variables')
xlabel('Time')
