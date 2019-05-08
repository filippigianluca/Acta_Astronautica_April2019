% plot belief 3D
% clear all; close all; clc

flag_input_list = 1;



%% u max mass
load('2sample1000');
LIST_F1 = LIST.table;
% stairs(LIST.F_Bel, LIST.Bel)
load('DataVolume_2samples1000');
% figure
% stairs(LIST.F_Bel, LIST.Bel)
LIST_F2 = LIST.table;




for i=1:length(LIST.table)-1
    max_FE_fun1(i) = LIST_F1{i+1,4};
    max_FE_fun2(i) = -LIST_F2{i+1,4};
    
    bpa_FE_fun1(i) = LIST_F1{i+1,6};
    bpa_FE_fun2(i) = LIST_F2{i+1,6};
end

x_min=12.2;
x_max=13.6;
y_min=600;
y_max=650;

N_nu=10000;

x = [12.3 un(max_FE_fun1,1e-6) 13.5];  
y = [599 uniquetol(max_FE_fun2,1e-8) 647];
[X,Y] = meshgrid(x, y);
F_bpa = zeros(size(X,1),size(X,2));


for j=1:length(x)
    for k=1:length(y)
        xx_yy = []; 
        xx_pos = find(max_FE_fun1(max_FE_fun1<x(j)));
        yy_pos = find(max_FE_fun2(max_FE_fun2>y(k)));
        xx_yy = intersect(xx_pos,yy_pos);
        if isempty(xx_yy)
            F_bpa(k,j) = 0;
        else
            F_bpa(k,j) = sum(bpa_FE_fun1(xx_yy));
        end
        
    end
end



contourf(X,Y,F_bpa,30)