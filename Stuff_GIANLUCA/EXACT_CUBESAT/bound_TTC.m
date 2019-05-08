function [lb_d, ub_d, lb_u, ub_u] = bound_TTC()


%% space_ttc: Communications system model
%
%   [M,P,info] = space_ttc(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = frequency, GHz
%       * x(2) = modulation,  
%                   PSK  = [0 1/8)
%                   BPSK = [1/8 2/8) 
%                   CFSK = [2/8 3/8) 
%                   BFSK = [3/8 4/8)
%                   FSK  = [4/8 5/8)
%                   DPSK = [5/8 6/8)
%                   QPSK = [6/8 7/8)
%                   NRZ  = [7/8 1]
%       * x(3) = antenna efficiency
%       * x(4) = antenna gain, dB
%       * x(5) = onboard loss, dB
%       * x(6) = other unmodelled losses, dB (polarization,implementation...)
%       * x(7) = mass of distribution network, kg
%       * x(8) = modulation index, rad
%       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
%
% * ep: Environmental parameters
%       * ep(1) = Bit Error Rate
%       * ep(2) = data volume, bits
%       * ep(3) = ground station G/T, dB
%       * ep(4) = range, km
%       * ep(5) = elevation angle, deg
%       * ep(6) = pointing accuracy, deg
%       * ep(7) = ground station antenna diameter, m
%       * ep(8) = ground station access time, min
%       * ep(9) = link margin, dB



lb_d = [ 7, ...   %       * x(1) = frequency, GHz
         0, ...   %       * x(2) = modulation,  
       0.4, ...   %       * x(3) = antenna efficiency
         1, ...   %       * x(4) = antenna gain, dB
        0.1, ...  %       * x(5) = onboard loss, dB
        0.5, ...  %       * x(6) = other unmodelled losses, dB (polarization,implementation...)
        0.1, ...  %       * x(7) = mass of distribution network, kg
        0.5, ...  %       * x(8) = modulation index, rad
          0, ...  %       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
         ];
     
     
ub_d = [ 10, ...   %       * x(1) = frequency, GHz
          1, ...   %       * x(2) = modulation,  
        0.6, ...   %       * x(3) = antenna efficiency
          5, ...   %       * x(4) = antenna gain, dB
          1, ...   %       * x(5) = onboard loss, dB
          2, ...   %       * x(6) = other unmodelled losses, dB (polarization,implementation...)
        0.5, ...   %       * x(7) = mass of distribution network, kg
          1, ...   %       * x(8) = modulation index, rad
          0.4         , ...   %       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
         ];

     
lb_u = [1e-6, ...       %       * ep(1) = Bit Error Rate
         1e5, ...       %       * ep(2) = data volume, bits
          20, ...       %       * ep(3) = ground station G/T, dB
         600, ...       %       * ep(4) = range, km
           5, ...       %       * ep(5) = elevation angle, deg
           5, ...       %       * ep(6) = pointing accuracy, deg
          10, ...       %       * ep(7) = ground station antenna diameter, m
           5, ...       %       * ep(8) = ground station access time, min
           5, ...       %       * ep(9) = link margin, dB
         ];

     
ub_u = [1e-4, ...       %       * ep(1) = Bit Error Rate
         1e7, ...       %       * ep(2) = data volume, bits
          40, ...       %       * ep(3) = ground station G/T, dB
         700, ...       %       * ep(4) = range, km
          15, ...       %       * ep(5) = elevation angle, deg
          15, ...       %       * ep(6) = pointing accuracy, deg
          20, ...       %       * ep(7) = ground station antenna diameter, m
          15, ...       %       * ep(8) = ground station access time, min
          15, ...       %       * ep(9) = link margin, dB
         ];
     
     
return