load('ro_sstl_200000feval_problem2.mat');
load('SSTL_bel_exact.mat')

cc = jet;
color_id_mu = 1+ceil(x(:,1)*63);
plot_belief(theta_list(:,end-1),theta_list(:,end))

figure(1)
hold on

for i = 1:size(fval,1)
    pf = plot(fval(i,1),1-fval(i,2),'*','Color',cc(color_id_mu(i),:)');
end
grid on;
xlabel('f_1 = cost(\mu,A)');
ylabel('1-f_2 = Bel(Pgen \geq Preq)');

c1=colorbar;
%.Label.string = '1-Bel(Pgen \geq Preq)'
c1.Label.String = '\mu';

% figure(2)
% hold on
% for i = 1:size(fval,1)
%     plot(fval(i,1),fval(i,2),'.','Color',cc(color_id_bel(i),:)');
%     
% end
% grid on;
% xlabel('f_1 = cost(\mu,A)');
% ylabel('f_2 = A');
% 
% c2=colorbar;
% c2.Label.String = 'Pl(Pgen \leq Preq)';