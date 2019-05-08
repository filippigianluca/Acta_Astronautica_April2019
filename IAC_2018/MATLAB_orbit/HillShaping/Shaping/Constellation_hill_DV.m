function Objective = Constellation_hill_DV(lambda,departure_elements,arrival_elements,TOF,mu,shape_flag,param_flag,du_i_h,du_others,flag)
% l = [a2 a4 a6 G2 G3 H2 H3 h2 h3]

lambda0 = [lambda(1) lambda(2) lambda(3) lambda(4) lambda(5)] ;
% tic
[c,A,b] = Shaping_functions(departure_elements,arrival_elements,shape_flag,param_flag,lambda0);
x  = A\(b-c);
[p] = Hill_Parameters_Ordering (x,lambda0,shape_flag,param_flag);


% t1 = tic;
[DeltaV,Time] = BinormalControl(departure_elements,arrival_elements,lambda(end),p,mu,du_i_h,du_others);
% t = toc(t1)
if flag == 0
    Objective = DeltaV  + abs(TOF-Time/86164.1004)*0.04;
elseif flag == 1
    Objective = DeltaV ;
end
% keyboard
end