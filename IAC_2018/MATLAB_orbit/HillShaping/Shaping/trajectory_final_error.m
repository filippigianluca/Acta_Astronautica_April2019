function error = trajectory_final_error(dep,arr,lambda,mu,nr,angle_rotation,flag)

% flag: 0 = Departure & Arrival Elements not rotated
% flag: 1 = Departure & Arrival Elements rotated


if nargin < 6 || flag == 1

    departure_elements = dep ; 
    arrival_elements = arr;
    
    
elseif flag == 0
    
    departure_elements = dep;
    departure_elements(5) = departure_elements(5) + angle_rotation;
    arr = wrapTo2Pi(arr);
    arrival_elements = RF_forward_hill_rotation(arr,dep(6),angle_rotation);
        
end
       
while arrival_elements(4) -2*pi*nr < departure_elements(4)
    arrival_elements(4) = arrival_elements(4) + 2*pi;
end

while abs(arrival_elements(5)-departure_elements(5)) > pi
    
    if arrival_elements(5) > departure_elements(5)
        arrival_elements(5) = arrival_elements(5) - 2*pi ;
    else
        arrival_elements(5) = arrival_elements(5) + 2*pi ;
    end
    
end
while abs(arrival_elements(6)-departure_elements(6)) > pi
    
    if arrival_elements(6) > departure_elements(6)
        arrival_elements(6) = arrival_elements(6) - 2*pi ;
    else
        arrival_elements(6) = arrival_elements(6) + 2*pi ;
    end
    
end

[c,A,b] = Shaping_functions(departure_elements,arrival_elements,1,0,lambda(1:5));
x  = A\(b-c);
[p] = Hill_Parameters_Ordering (x,lambda(1:5),1,0);
[r,G] = Hill_parameters_analytical_shaping(p,departure_elements(4),1);
[his,F_n,Inclination,dF_n,dInclination,~] = BinormalControl_analytical(departure_elements,arrival_elements,r,G,lambda(6),p,5000);
[F_r,F_u] = Hill_control_analytical_shaping3(F_n,dF_n,Inclination,dInclination,mu,p,departure_elements(4),1) ;
U = u_constraction_analytical(his,1500);

[~,Hill_elements] = RK4_Hill (U,departure_elements,F_r,F_u,F_n,mu);
if size(Hill_elements,2) == 6
    Hill_elements = Hill_elements';
end
arrival_reconstructed_elements = Hill_elements(:,end);

if nargin < 6 || flag == 1

    error = ( arrival_reconstructed_elements - arr ) ./ arr *100 ;
    
elseif flag == 0
    
    arrival_reconstructed_elements = RF_backward_hill_rotation(arrival_reconstructed_elements,dep(6),angle_rotation);
    error = ( arrival_reconstructed_elements - arr ) ./ arr *100 ;
    
end





end