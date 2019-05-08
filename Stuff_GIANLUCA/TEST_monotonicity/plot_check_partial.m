clear all; close all; clc


color = {'r', 'k'};
N = [20000];% 30000];


a=1;
figure
for i = N
load(strcat('solution_max_ui_uij_fixed_',num2str(i),'feval'))

% curva 1
hold on

plot([1:4],F_MAX(2,1:4),color{a})
plot([1:4],Partial_curve{1, 2}.Belief_FE_function(2:end),'b')

scatter([1:4],F_MAX(2,1:4),'filled',color{a})
scatter([1:4],Partial_curve{1, 2}.Belief_FE_function(2:end),'filled','b')
a=a+1;

end
legend('max_{ui}','partial curves')


figure
for i = N
load(strcat('solution_max_ui_uij_fixed_',num2str(i),'feval'))

% curva 2
hold on

plot([1:8],F_MAX(3,1:8),'r')
plot([1:8],Partial_curve{1, 3}.Belief_FE_function(2:end),'b')

scatter([1:8],F_MAX(3,1:8),'filled','r')
scatter([1:8],Partial_curve{1, 3}.Belief_FE_function(2:end),'filled','b')
end
legend('max_{ui}','partial curves')

