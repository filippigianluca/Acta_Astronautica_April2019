function [F] = function_F(d,u,par)
   
%%[2 2 2]
f_1 = 10*u(1)^2 + abs(u(2))*u(5)^2 + u(6)^4/100 +d(1)^2*abs(d(2));        
f_2 = abs(u(3)) + u(4)^2*abs(u(5))/10 + u(6)^2 + abs(d(1));   




F = f_1 + f_2;                                                        


global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end