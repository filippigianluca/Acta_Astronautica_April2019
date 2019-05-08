function U = u_constraction(x,u,h_Fru,h_Fn)


U = [];
if x(1,1) > u(1)+1e-3
    
    n_steps = ceil((U(1,1)-h_Fru - u(1))/h_Fru +1) ;    
    U = [U ; linspace(u(1),U(1,1)-h_Fru,n_steps)' ];
        
end

l = size(x,1); 

for ii = 1:l
   keyboard
    n_steps = ceil((U(ii,2) -U(ii,1))/h_Fn +1) ;    
    U = [U ; linspace(U(ii,1),U(ii,2),n_steps)' ];
    
    if ii ~= l
        
            n_steps = ceil((U(ii+1,1) - U(ii,2)-2*h_Fru)/h_Fru +1) ;    
            U = [U ; linspace(U(ii,2)+h_Fru,U(ii+1,1)-h_Fru,n_steps)' ];
            
    else
        if U(ii,2) < u(2)
            n_steps = ceil((u(end) - U(ii,2)-h_Fru)/h_Fru +1) ;    
            U = [U ; linspace(U(ii,2)+h_Fru,u(end),n_steps)' ];
        end
    end       
       
    
    
end



end