

rip_tot = 10000;


[state_minmax,  s_minmax]   = histogram_minmax_2(rip_tot);

[state_nominal, s_nominal] = histogram_nominal_2(rip_tot);




%%
figure
hold on

h1_nominal = hist(state_nominal,[1 2]);
h1_nominal = h1_nominal/sum(h1_nominal);
% bar(h1, 'DisplayName', 'N transitions');

% legend('show');
% title('min-max')



%%
h1_minmax = hist(state_minmax,[1 2]);
h1_minmax = h1_minmax/sum(h1_minmax);
% bar(h2, 'DisplayName', 'Time');
% legend('show');
% title('min-max')


bar([1 2],[h1_nominal;h1_minmax]')
legend('show');

