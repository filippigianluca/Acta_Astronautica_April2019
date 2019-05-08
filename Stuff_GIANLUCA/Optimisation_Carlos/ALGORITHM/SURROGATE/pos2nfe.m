function [ nfe ] = pos2nfe( pos, n_int )

nfe = 0;
num_divisioni_vecchie = 1;

for i = 1:length(n_int)-1 

    num_divisioni_vecchie = num_divisioni_vecchie*n_int(i);

    nfe = nfe + (pos(i+1)-1)*num_divisioni_vecchie;     % indici: matrice --> vettore num_divisioni_vecchie(i)

end
nfe = nfe + pos(1); 


end

