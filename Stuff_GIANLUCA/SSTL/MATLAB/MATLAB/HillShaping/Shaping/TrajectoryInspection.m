function [U,x] = TrajectoryInspection (u0,uf,bound)

% x = 0 -> Line of node passage where i change
% x = 1 -> Normal to line of node passage where h change

nr = floor((uf-u0)/(2*pi));

r2d = 180/pi ;
d2r = pi/180 ;

aux = [floor(u0*r2d) : 1 : ceil(uf*r2d)]';
frac = (aux/90);
R = rem(aux,90);

index = find(R == 0);

x = mod(frac(index),2);

U = [aux(index)*d2r-bound , aux(index)*d2r+bound];

if u0 > U(1,1) && u0 < U(1,2)
    U (1,1) = u0 ; 
    
end

if uf > U(end,1) && uf < U(end,2)
%     U (end,2) = uf ;     
    U = U(1:end-1,:); % Non considerare ultimo pezzo perchE' cambia inclinazione e RAAN
    x = x(1:end-1);
end

if U(1,1)~=u0
    if (u0 < bound && u0 > 0)
        
        U = [u0 bound ; U];
        x = [0 ; x];
        
    elseif (u0 < pi/2+bound && u0 > pi/2)
        
        U = [u0 pi/2+bound ; U];
        x = [1 ; x];
        
    elseif (u0 < pi+bound && u0 > pi)
        
        U = [u0 pi+bound ; U];
        x = [0 ; x];
        
    elseif (u0 < 3/2*pi+bound && u0 > 3/2*pi)
        
        U = [u0 pi*3/2+bound ; U];
        x = [1 ; x];
        
    end
end

% uft = wrapTo2Pi(uf);
% 
% if U(end,2)~=uf
%     if (uft > pi/2-bound && uft < pi/2)
%         
% %         U = [U ; (pi/2+nr*2*pi)-bound uf];
% %         x = [x ; 1];
% 
%         
%     elseif (uft > pi-bound && uft < pi)
%         
% %         U = [U ; (pi+nr*2*pi)-bound uf];
% %         x = [x ; 0];
% 
%         
%     elseif (uft > 3/2*pi-bound && uft < 3/2*pi)
%         
% %         U = [U ; (3/2*pi+nr*2*pi)-bound uf];
% %         x = [x ; 1];
% 
%         
%     elseif (uft > 2*pi-bound && uft < 2*pi)
%         
% %         U = [U ; (2*pi+nr*2*pi)-bound uf];
% %         x = [x ; 0];
% 
%     end
% end



end