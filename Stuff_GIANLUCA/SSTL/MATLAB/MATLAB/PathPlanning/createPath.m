function [r, dr, ddr, timespan, a] = createPath(r0, v0, rf, vf, af, t0, tf, dt)
% Create Path Planning based on a 7th polynomial
% 
% NOTE: 
%  It is assumed that intial t0 is zero for now. 
%
% INPUTS:
%  
%     r0  - Initial Position (size: 3) 
%     v0  - Initial Velocity (size: 3) 
%     rf  - Final Position (size: 3) 
%     vf  - Final Velocity (size: 3) 
%     af  - Final Acceleration (size: 3) 
%     t0  - Initial Time (IT IS ASSUMED THAT t0 IS ZERO)
%     tf  - Final Time 
%     dt  - Time Step 
%
% OUTPUTS:
%      r  - Position Profile (size: 3 x n) 
%      v  - Velocity Profile (size: 3 x n)
%      a  - Acceleration Profile (size: 3 x n)
%
% AUTHOR:
%  (C) Juan Manuel Romero Martin -- 2016
%==========================================================================


%==========================================================================
% STEP 1: PERFORM INPUTS SANITY CHECKS 
%==========================================================================

if( length(r0) ~= length(v0) )
    warning('The r0 and v0 have to be with the same dimension');
end

if( length(rf) ~= length(vf) )
    warning('The rf and vf have to be with the same dimension');
end

if( length(r0) ~= length(rf) )
    warning('The r0 and rf have to be with the same dimension');
end

if( length(v0) ~= length(vf) )
    warning('The v0 and vf have to be with the same dimension');
end

if( length(af) ~= length(vf) )
    warning('The af and vf have to be with the same dimension');
end

%==========================================================================
% STEP 2: SOLVE THE POLYNOMIAL PATH 
%==========================================================================

% Create the partial Matrix for the provided boundaries conditions
A_part = [1   0      0       0         0;  ...% Element i for r0  
          0   1      0       0         0;  ...% Element i for v0
          1  tf   tf^2    tf^3      tf^4;  ...% Element i for rf
          0   1   2*tf  3*tf^2    4*tf^3;  ...% Element i for vf
          0   0      2    6*tf   12*tf^2];    % Element i for af
     
% Solve for each of the dimenstions
for i = 1 : length(r0)
    b = [r0(i) v0(i) rf(i) vf(i) af(i)]';
    a_n5(:,i) = A_part\b; 
end

% Fill the rest of the matrix to meet the order of the polynomial
% @TODO: Modify this to be dynamic to the degree of the polynomial
% A = [a_n5; zeros(size(t_vec,1)-5, 3)]; 
a = [a_n5; zeros(3, 3)]; 

%==========================================================================
% STEP 3: COMPUTE ALL THE COMPONENTS FOR THE WHOLE TIME-SPAN
%========================================================================== 

% Compute the Number of points for the time-span
num_points = floor((tf - t0)/dt);

% Compute the Time-span
timespan  = linspace(t0, tf, num_points);
t = timespan; % to simplify the following statements

% WARNING: t-vectors must have n rows for n-order polynomials.
t_vec    = [ ones(1,num_points);                   t;                  t.^2;    t.^3;    t.^4;    t.^5;    t.^6;    t.^7];  % For Position
t_vec_p  = [zeros(1,num_points);  ones(1,num_points);                2*t.^1;  3*t.^2;  4*t.^3;  5*t.^4;  6*t.^5;  7*t.^6];  % For Velocity
t_vec_pp = [zeros(1,num_points); zeros(1,num_points);  2*ones(1,num_points);  6*t.^1; 12*t.^2; 20*t.^3; 30*t.^4; 42*t.^5];  % For Acceleration

% Compute the Trajectory
r   = a'*t_vec; 
dr  = a'*t_vec_p;
ddr = a'*t_vec_pp;

end % End Function