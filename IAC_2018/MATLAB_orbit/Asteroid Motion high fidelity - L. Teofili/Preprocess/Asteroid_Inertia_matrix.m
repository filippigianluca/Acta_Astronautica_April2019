function [I,volume,mass]=Asteroid_Inertia_matrix(r1,r2,r3,rho)
%Asteroid_Inertia_matrix computes principal inertia matrix, the volume and the mass 
%  of the Asteroid. 
%  [I,volume,mass]=Asteroid_Inertia_matrix(r1,r2,r3,rho) as input the vector
%  vertex matrices r1 r2 r3 [given by load_model function] and the mean 
%  asteroid density and give back principal inertia matrix, the volume and the mass 
%  of the Asteroid. 
%  References are "Explicit exact formulas for the 3-d tetrahedron inertia 
%  tensor in terms of its vertex coordinates", F. Tonon.
%  Lorenzo Teofili, University of Rome La Sapienza 
%  Nicolas Thiry, University of Strathclyde, 5/11/2015

[l,m] = size(r1);

volume = 0;

I=zeros(3,3);

for i=1:l
    
    %Inertia computation
    x1 = r1(i,1);
    x2 = r2(i,1);
    x3 = r3(i,1);
    y1 = r1(i,2);
    y2 = r2(i,2);
    y3 = r3(i,2);
    z1 = r1(i,3);
    z2 = r2(i,3);
    z3 = r3(i,3);
    
    J  = (det([x1 x2 x3;y1 y2 y3;z1 z2 z3]));                                                           % WTF?? Is that the Volume? if we need a tetraedron we have to divide per 6!
    a  = (y1^2 + y1*y2 + y2^2 + y1*y3 + y2*y3 + y3^2 + z1^2 + z1*z2 + z2^2 + z1*z3 + z2*z3 + z3^2)/60;  % Why 60????? Why the mixed product?
    b  = (x1^2 + x1*x2 + x2^2 + x1*x3 + x2*x3 + x3^2 + z1^2 + z1*z2 + z2^2 + z1*z3 + z2*z3 + z3^2)/60;
    c  = (x1^2 + x1*x2 + x2^2 + x1*x3 + x2*x3 + x3^2 + y1^2 + y1*y2 + y2^2 + y1*y3 + y2*y3 + y3^2)/60;
    a_ = (2*y1*z1 + y2*z1 + y3*z1 + y1*z2 + 2*y2*z2 + y3*z2 + y1*z3 + y2*z3 + 2*y3*z3)/120;
    b_ = (2*x1*z1 + x2*z1 + x3*z1 + x1*z2 + 2*x2*z2 + x3*z2 + x1*z3 + x2*z3 + 2*x3*z3)/120;
    c_ = (2*x1*y1 + x2*y1 + x3*y1 + x1*y2 + 2*x2*y2 + x3*y2 + x1*y3 + x2*y3 + 2*x3*y3)/120;
    
    I = I + rho*J*[a  -b_ -c_ ; ... 
                  -b_  b  -a_ ; ...
                  -c_ -a_  c ];
              
    volume = volume + J'/6; % also J/6
              
end

mass = volume*rho;

end
    