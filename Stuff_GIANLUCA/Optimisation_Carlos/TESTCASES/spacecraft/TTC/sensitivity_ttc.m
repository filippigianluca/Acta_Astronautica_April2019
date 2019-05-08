% sensitivity analysis of space_ttc

clear all;close all;clc;
x=rand(1,9);
% x(2)=rand;
ep = ones(1,9);
[M,P,info] = space_ttc(x,ep);

% for x
for i=1:9
idx=1;   
    for x_value=1:20
        x(i) =x_value;
        if i ==2
            x(i)=x_value/20;
        end
        [M,P,info] = space_ttc(x,ep);
        Pot(idx)=P;
        id(idx)=x_value;
        idx=idx+1;
    end
    figure
    plot(id(:),Pot)
    hold off
    legend(strcat('x', num2str(i)))
end

% % for ep
% for i=1:9
%     
%     for idx=1:20
%         ep(i) =idx;
%         [M,P,info] = space_ttc(x,ep);
%         Pot(idx)=P;
%     end
%     figure
%     plot([1:20],Pot)
%     hold off
%     legend(strcat('ep', num2str(i)))
% end