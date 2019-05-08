% close all; clear all; clc;

%% INPUT values
num_punti_doe = 10;
dim_u = 1;                            
dim_d = 1;
n_obj = 2;

lb_u = {repmat({[-5,-3,-1]},[dim_u,1]),repmat({[-5,-3,-1]},[dim_u,1])};
ub_u = {repmat({[-4,0,3]},[dim_u,1]),repmat({[-4,0,3]},[dim_u,1])};

%% DOE generation
x_doe_aux = rand(num_punti_doe,dim_u+dim_d);
f_doe_aux = 10*rand(num_punti_doe,1);


% Surrogate structure generation

surrogate = cell(1,n_obj);
for obj = 1:n_obj
    map_u_info = get_map_info(lb_u{obj}, ub_u{obj});
    
    num_fe = 1;
    for dim = 1:dim_u
        num_fe = num_fe * map_u_info.n_int{dim};
    end

    surrogate{obj} = struct;
    surrogate{obj}.dim_d = dim_d;
    surrogate{obj}.dim_u = dim_u;
    surrogate{obj}.method = 'kriging';
    surrogate{obj}.corrfun = @corrgauss;
    surrogate{obj}.regrfun = @regpoly0;
    surrogate{obj}.training = str2func([lower(surrogate{obj}.method) '_training']);
    surrogate{obj}.predictor = str2func([lower(surrogate{obj}.method) '_predictor']);
    surrogate{obj}.indicator = str2func([lower(surrogate{obj}.method) '_EI']);
    surrogate{obj}.model = cell(1,num_fe);
    surrogate{obj}.x_doe = cell(1,num_fe);
    surrogate{obj}.f_doe = cell(1,num_fe);
    surrogate{obj}.map_info = map_u_info;
end


%% FUNCTION
for obj = 1:n_obj
    [surrogate{obj}] = surrogate_add(x_doe_aux,f_doe_aux,surrogate{obj});
end