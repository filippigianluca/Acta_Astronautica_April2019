function [cost] = sstl2c(d);

    area1 = 0.0007996; %0.0007996; !!!
    area2 = 0.0026639; %0.0007996; !!!p
    C1 = 70/area1;
    C2 = 350/area2; 
    A = d(2);
    mu = d(1);
    cost = A*(mu*C1+(1-mu)*C2);

return