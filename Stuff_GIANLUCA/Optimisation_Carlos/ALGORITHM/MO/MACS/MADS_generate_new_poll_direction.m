function b = MADS_generate_new_poll_direction(l,n)

% following exactly Mesh Adaptive Direct Search Algorithms for Constrained
% Optimisation

ihat = randi(n);
set = -2^l+1:2^l-1;
b = zeros(n,1);
plus_or_min = randi(2);
b(ihat) = (-1)^plus_or_min*2^l;

for i = 1:ihat-1
    
    b(i) = set(randi(length(set)));
    
end

for i = ihat+1:n
    
    b(i) = set(randi(length(set)));
    
end

end