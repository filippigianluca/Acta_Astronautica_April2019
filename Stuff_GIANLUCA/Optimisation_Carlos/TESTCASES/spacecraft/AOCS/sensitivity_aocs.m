% sensitivity analysis of space_aocs

clear all;close all;clc;
x=rand(1,14);
ep = rand(1,13);
idx=1;

% % for x
for i=1:14
    
    for idx=1:20
        x(i) =idx;
        [M,P,info] = space_aocs(x,ep);
        Pot(idx)=P;
    end
    figure
    plot([1:20],Pot)
    hold off
    legend(strcat('x', num2str(i)))
end

% x_2 = 0.005:0.0005:0.02;
% x_3 = 0.034:0.01:0.0885;
% x_4 = 0.034:0.01:0.0885;
% x_5 = 0.5:0.05:0.7;
% 
% for index =x_2;
%         hold on
%         x(2) = index;
%         [M,P,info] = space_aocs(x,ep);
%         Pot(idx)=P;
%         idx = idx+1;
% end
% figure
% plot(x_2,Pot)
% hold off
% legend(strcat('x', num2str(i)))


% for ep
% for i=1:13
%     
%     for idx=1:20
%         ep(i) =idx;
%         [M,P,info] = space_aocs(x,ep);
%         Pot(idx)=P;
%     end
%     figure
%     plot([1:20],Pot)
%     hold off
%     legend(strcat('ep', num2str(i)))
% end