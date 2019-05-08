% plot for different thresholds in the constraint the minmax solution: 
% 1. d-minmax, u-minmax
% 2. d-minmax, u-max-constraint

%% defin settings

nu_vec = [500 600 700 800];
func_handle = @Acta_Astronautica_objconstr1;
type_plot = 2; % 1. two points
               % 2. one point

%% run

lineStyles=linspecer(length(nu_vec),1); 
par.fix.time = 365;
hold on


index_plot = 1;
for nu = nu_vec
    plot([nu nu], [10 18],'k','linewidth',1);
    par.fix.nu = nu;
    load(strcat('Acta_nu',num2str(nu),'_so_nfeval_minmax_2000000nfeval_inner_outer_15000'));
    
    f_minmax_archive = cell2mat(minmax.archiveG_cc(2:end,4));
    f_minmax = minmax.f{1};
    pos = find(f_minmax_archive==f_minmax);
    d_minmax = minmax.d;
    u_minmax = minmax.u{1};
    u_max_c  = minmax.archiveG_cc{pos+1,3};
    
    load('map_u_info');
    u_max_c = map_affine(u_max_c, map_u_info);
    
    c_max    = nu - minmax.archiveG_cc{pos+1,6};
    
    output = func_handle(d_minmax,u_max_c,par);
    mass_c_max = output.f;
    
    output = func_handle(d_minmax,u_minmax,par);
    c_minmax = nu - output.c;
    
    if type_plot == 1
    
        s_max_c(index_plot) = scatter(c_max, mass_c_max,'filled','DisplayName',strcat('constr SO max violation'));
        s_max_c(index_plot).MarkerFaceColor = lineStyles(index_plot,:);
        s_max_c(index_plot).SizeData = 300;
        s_max_c(index_plot).Marker = 'd';
        
        s_minmax(index_plot) = scatter(c_minmax, f_minmax,'filled','DisplayName',strcat('constr SO'));
        s_minmax(index_plot).MarkerFaceColor = lineStyles(index_plot,:);
        s_minmax(index_plot).SizeData = 300;    
        
        p(index_plot) = plot([c_max c_minmax], [mass_c_max f_minmax],'Color',lineStyles(index_plot,:),'linewidth',1,'DisplayName',strcat('constr min-max \nu = ',num2str(nu)));
    
    
    else
        
        s_max_c(index_plot) = scatter(c_max, f_minmax,'filled','DisplayName',strcat('constr min-max \nu = ',num2str(nu)));
        s_max_c(index_plot).MarkerFaceColor = lineStyles(index_plot,:);
        s_max_c(index_plot).SizeData = 300;     
        
    end
    
    
    index_plot = index_plot+1;
end

grid on
xlabel('f_V (GB)')
ylabel('Mass (kg)')

title(strcat('constrained min-max pareto front'))

if type_plot == 1
    lgd = legend(p(1,:));  
else
    lgd = legend(s_max_c(1,:));  
end
