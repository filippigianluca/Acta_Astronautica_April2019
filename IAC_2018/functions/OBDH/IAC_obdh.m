function [M, P, info] =  IAC_obdh(x, ep)

% (1) http://propagation.ece.gatech.edu/ECE6390/project/Sum2015/team5/satellite-downlink.html
%
% compression to JPEG: high quality (10:1) for land (29.2%) and low quality
% (50:1) for water (70.8%). 



% OUTPUT: - M, mass of OBDH
%         - P, power of OBDH
%         - info, - data volume after compression



%% ------------------------------------------------------------------------



% design
TIME = x(1);  % variable

if x(2) < 0.25 
    obdh = 'vectronic'; % https://www.vectronic-aerospace.com/space-applications/payload-data-handling-system-vpdhs/
elseif x(2) < 0.5
    obdh = 'invented1';
elseif x(2) < 0.75
    obdh = 'invented2';
else
    obdh = 'invented3';
end



% epistemic
N_foto_to_store  = ep(1); % maximum number of pictures to compress and store in all he loops
N_foto_tot       = ep(2); 
V_to_store       = ep(3);
V_tot            = ep(4);
M_margin       = (1+ep(5)/100);
P_margin       = (1+ep(6)/100);

    






%% ------------------------------------------------------------------------

% image compression to JPEG from the link (1)
WR = 0.708; % water_rate
CW = 0.02;  % compression_rate_water
LR = 0.292; % land_rate
CL = 0.1;   % compression_rate_land

C = WR*CW+LR*CL; % Total Compressed Daily Data Size

% data volume after compression TO STORE
V_compressed_to_storage = V_to_store*C;

% data volume after compression TOTAL
V_compressed_tot = V_tot*C;



%
% 1	Urban	1%	10:1
% 2	Land/non-urban	28%	20:1
% 3	Water	71%	50:1

switch obdh
    case 'vectronic' 
    M = 2.3;          % (kg)
    P = 15;           % (Watt)
    info.storage = 4; % (Gbytes)
    
    case 'invented1'
    M = 2;          % (kg)
    P = 20;           % (Watt)
    info.storage = 4; % (Gbytes)   
    
    case 'invented2'
    M = 1.5;          % (kg)
    P = 22;           % (Watt)
    info.storage = 4; % (Gbytes)    
    
    case 'invented3'
    M = 3;            % (kg)
    P = 30;           % (Watt)
    info.storage = 4; % (Gbytes)   
end

N_obdh = 1;
storage_start = info.storage;
while V_compressed_to_storage > info.storage
    N_obdh = N_obdh + 1;
    info.storage = info.storage + storage_start;
end
M = M*N_obdh;
P = P*N_obdh;



% power increase due to: 
% 1. compression
% 2. storage
P = P + (0.000001*N_foto_tot + 0.000001*N_foto_to_store)*P_margin;
M = M + (0.0000001*N_foto_tot + 0.000001*N_foto_to_store)*M_margin;





%--------------------------------------------------------------------------
info.V_compressed_tot = V_compressed_tot;


end