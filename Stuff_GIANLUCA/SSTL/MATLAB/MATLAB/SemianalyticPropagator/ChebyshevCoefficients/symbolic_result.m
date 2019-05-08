function P_final = symbolic_result(c,n,a,b)

% Marilena Di Carlo, 2015

if a == -1  && b == 1
    x = sym('x');
else
    y = sym('y');
    x = -1 + 2 * y /(b-a) - 2 * a / (b-a);
end

T = [1; x; zeros(n-1,1)];

for i = 2 : n
    
 T(i+1) = 2*x*T(i) - T(i-1);
 
end

% keyboard
T = collect(T);

P_final = sum(c * T);


end