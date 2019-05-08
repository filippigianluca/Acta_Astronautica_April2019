% plot results of minmax problem:
clear all; close all; clc

load('IAC2018_nfeval_1e6nInnerOuter_2e4_npop_3sqp_end')
% load('IAC2018_1e6_2e3_2e3_end');


init_tc   = str2func(strcat('init_tc_mo_IAC2018'));
[ problem_minmax ] = init_tc();


for i=2:size(minmax.output.ARCHIVE{1, 1},1)
    f_minmax(i-1) = minmax.output.ARCHIVE{1, 1}{i,4};
    N_minmax(i-1) = minmax.output.ARCHIVE{1, 1}{i,7};
    
    d_minmax = minmax.output.ARCHIVE{1, 1}{i,1};
    u_minmax = minmax.output.ARCHIVE{1, 1}{i,2};
    MASS(i-1) = CUBESAT_5subsystems_MASS(d_minmax, u_minmax, problem_minmax);
    RES(i-1) = - CUBESAT_5subsystems_RES(d_minmax, u_minmax, problem_minmax);
end

figure
subplot(3,1,1);
scatter(N_minmax, f_minmax, 'filled','b')
hold on
plot(N_minmax, f_minmax,'b')
xlabel('nfeval')
ylabel('F=M/RES')

subplot(3,1,2);
scatter(N_minmax(2:end), MASS(2:end), 'filled','b')
hold on
plot(N_minmax(2:end), MASS(2:end),'b')
xlabel('nfeval')
ylabel('mass')

subplot(3,1,3);
scatter(N_minmax(2:end), RES(2:end), 'filled','b')
hold on
plot(N_minmax(2:end), RES(2:end),'b')
xlabel('nfeval')
ylabel('RES')

