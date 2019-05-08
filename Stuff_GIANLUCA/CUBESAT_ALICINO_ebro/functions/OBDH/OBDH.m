function R = OBDH(d, u, par)


% OUTPUT: data rate R




TIME    % variable
Tac     % from payload
T_orbit % from payload
N_foto  % from payload
V_bits_per_foto 




V_bits = V_bits_per_foto*N_foto;

% data rate R: V is the data volume, in bits, to be transmitted, and Tac is the access time, in seconds, to the ground station
R = V_bits/Tac; 
end