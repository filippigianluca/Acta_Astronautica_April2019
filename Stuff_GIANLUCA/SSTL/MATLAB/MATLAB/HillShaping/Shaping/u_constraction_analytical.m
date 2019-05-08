function u = u_constraction_analytical(x,n)

l = length(x); 
u = [];
for ii = 1:l-1
   
    if x(ii+1)<x(ii)+0.1
        u = [u ; x(ii)];
    else
        aux = linspace(x(ii),x(ii+1),n)';
        u = [u ; aux(1:end-1)] ;
    end
    
end

u = [u ; x(end)] ;


end