% plot belief 3D
% clear all; close all; clc

flag_input_list = 1;
flag_method = 2;



%% u max mass
load('2sample1000');
LIST_F1 = LIST.table;
load('DataVolume_2samples1000');
LIST_F2 = LIST.table;




for i=1:length(LIST.table)-1
    max_FE_fun1(i) = LIST_F1{i+1,4};
    max_FE_fun2(i) = -LIST_F2{i+1,4};
    
    bpa_FE_fun1(i) = LIST_F1{i+1,6};
    bpa_FE_fun2(i) = LIST_F2{i+1,6};
end




[max_fun1_sort, n_fun1_sort] = sort(max_FE_fun1);
[max_fun2_sort, n_fun2_sort] = sort(max_FE_fun2);






% %-------------
% max_fun2_sort_opp = max_fun2_sort(end:-1:1);
% [X_opp,Y_opp] = meshgrid(max_fun1_sort, max_fun2_sort);
% %-------------------------


if flag_method == 1
    max_fun1_sort = [12.5 : 0.01 : max_fun1_sort(1)-0.1 max_fun1_sort max_fun1_sort(end)+0.1:0.01:13.5];
    max_fun2_sort = [390 : max_fun2_sort(1)-1 max_fun2_sort max_fun2_sort(end)+1:445];
    
    [X,Y] = meshgrid(max_fun1_sort, max_fun2_sort);
    
    F = zeros(size(X,1),size(X,2));
    for j=1:size(X,1)
        for k=1:size(X,2)
            
            for m = 1: length(max_FE_fun1)
                
                if max_FE_fun1(m) <= X(j,k) && max_FE_fun2(m) >= Y(j,k)
                    F(j,k) = F(j,k) + bpa_FE_fun1(m);
                end
            end
        end
    end
    
    
    
    
elseif flag_method == 2
    
    [X,Y] = meshgrid(max_fun1_sort, max_fun2_sort);
    
    F = zeros(size(X,1),size(X,2));
    
    
 for ii = 1:length(n_fun1_sort)
    
    pos_FEs = n_fun1_sort(1:ii);
    bpa_fun2_sort = zeros(1,length(n_fun1_sort));
    bpa_fun2_sort(pos_FEs) = bpa_FE_fun2(n_fun1_sort(pos_FEs));
    
    a = cumsum(bpa_fun2_sort);
    F(:,ii) = a(end:-1:1);




 end
    
 
 
 elseif flag_method == 3
     position = 1:length(max_fun1_sort);
%      position = 100:100:length(max_fun1_sort);
     max_fun1_sort_new = max_fun1_sort(position);     
     max_fun2_sort_new = max_fun2_sort(position);   
    if  max_fun1_sort(end) ~= max_fun1_sort_new(end)
        max_fun1_sort_new = [max_fun1_sort_new max_fun1_sort(end)];
        max_fun2_sort_new = [max_fun2_sort_new max_fun1_sort(end)];
        position = [position length(max_fun1_sort)];
    end
     
    n_fun1_sort = n_fun1_sort(position);
    [~,n_fun1_sort] = sort(n_fun1_sort);
     bpa_FE_fun1_new(1) = sum(bpa_FE_fun1(1:position));
     bpa_FE_fun2_new(1) = sum(bpa_FE_fun2(1:position));

     for iii =2:length(max_fun1_sort_new)
         bpa_FE_fun1_new(iii) = sum(bpa_FE_fun1(1+position(iii-1):position(iii)));
         bpa_FE_fun2_new(iii) = sum(bpa_FE_fun2(1+position(iii-1):position(iii)));

     end
     
    max_fun1_sort = max_fun1_sort_new;
    max_fun2_sort = max_fun2_sort_new;
    bpa_FE_fun1 = bpa_FE_fun1_new;
    bpa_FE_fun2 = bpa_FE_fun2_new;

    
    [X,Y] = meshgrid(max_fun1_sort, max_fun2_sort);
    
    F = zeros(size(X,1),size(X,2));
    
    
 for ii = 1:length(n_fun1_sort)
    
    pos_FEs = n_fun1_sort(ii);
    bpa_fun2_sort = zeros(1,length(n_fun1_sort));
    bpa_fun2_sort(pos_FEs) = bpa_FE_fun2(n_fun1_sort(pos_FEs));
    
    a = cumsum(bpa_fun2_sort);
    if ii >1
        F(:,ii) = F(:,ii-1)' + a(end:-1:1);
    else
        F(:,ii) = a(end:-1:1);
    end
 end

end



contourf(X,Y,F)