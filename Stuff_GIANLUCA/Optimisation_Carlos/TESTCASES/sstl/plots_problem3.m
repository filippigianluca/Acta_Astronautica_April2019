% load('ro_sstl_200000feval_problem3.mat');
load('ro_sstl_p3orig_250000feval.mat');

cc = jet;
color_id_bel = 1+ceil(fval(:,3)*63);
color_id_mu = 1+ceil(x(:,1)*63);
figure(1)
hold on
for i = 1:size(fval,1)
    plot3(fval(i,1),fval(i,2),fval(i,3),'.','Color',cc(color_id_mu(i),:)');
    
end
grid on;
xlabel('f_1 = cost(\mu,A)');
ylabel('f_2 = A');
zlabel('Pl(Pgen \leq Preq)');

c1=colorbar;
%.Label.string = '1-Bel(Pgen \geq Preq)'
c1.Label.String = '\mu';

figure(2)
hold on
for i = 1:size(fval,1)
    plot(fval(i,1),fval(i,2),'.','Color',cc(color_id_bel(i),:)');
    
end
grid on;
xlabel('f_1 = cost(\mu,A)');
ylabel('f_2 = A');

c2=colorbar;
c2.Label.String = 'Pl(Pgen \leq Preq)';