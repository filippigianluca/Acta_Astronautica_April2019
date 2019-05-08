function [lb_d, ub_d, lb_u, ub_u] = bound_TTC_Alicino_Paper()


%% space_ttc: Communications system model
%
%   [M,P,info] = space_ttc(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = frequency, GHz                           (ALICINO design       [7 10])
%       * x(2) = modulation,                              (ALICINO design       [0 1])
%       * x(3) = antenna efficiency                       (ALICINO uncertainty  [0.6 0.8] [0.8 0.9])
%       * x(4) = antenna gain, dB                         (ALICINO uncertainty  [1 3][3 5] trasmit AG)
%       * x(5) = onboard loss, dB                         (ALICINO uncertainty  [0.1 0.5][0.5 1])
%       * x(6) = other unmodelled losses, dB              (ALICINO uncertainty  [0.5 1.5][1.5 2])
%       * x(7) = mass of distribution network, kg         (ALICINO uncertainty  [0.1 0.3][0.2 0.5])
%       * x(8) = modulation index, rad
%       * x(9) = amplifier type,                          (ALICINO design       [0 1])
%                   
%
% * ep: Environmental parameters
%       * ep(1) = Bit Error Rate                          (ALICINO fix, 1e-5)
%       * ep(2) = data volume, bits                       (ALICINO fix, 1e6)
%       * ep(3) = ground station G/T, dB                  (ALICINO fix, 30)
%       * ep(4) = range, km                               (ALICINO fix, 640)
%       * ep(5) = elevation angle, deg                    (ALICINO fix, 10)
%       * ep(6) = pointing accuracy, deg                  (ALICINO fix, 5)
%       * ep(7) = ground station antenna diameter, m      (lambda/pi * (10e6/0.55)^0.5; lambda = c/(1e9*f))
%       * ep(8) = ground station access time, min         (ALICINO fix, 10)
%       * ep(9) = link margin, dB                         (margin on the P, i can fix to 0)



lb_d = [ 7, ...   %       * x(1) = frequency, GHz                                                (ALICINO design)
         0, ...   %       * x(2) = modulation,                                                   (ALICINO design)
       0.4, ...   %       * x(3) = antenna efficiency                                            (ALICINO uncertainty)                   
         1, ...   %       * x(4) = antenna gain, dB                                              (ALICINO uncertainty)
        0.1, ...  %       * x(5) = onboard loss, dB                                              (ALICINO uncertainty)
        0.5, ...  %       * x(6) = other unmodelled losses, dB (polarization,implementation...)  (ALICINO uncertainty)
        0.1, ...  %       * x(7) = mass of distribution network, kg                              (ALICINO uncertainty)
        0.5, ...  %       * x(8) = modulation index, rad
          0, ...  %       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]                (ALICINO design)
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

     
lb_u = [1e-6, ...       %       * ep(1) = Bit Error Rate                                          (ALICINO fix, 1e-5)
         1e5, ...       %       * ep(2) = data volume, bits                                       (ALICINO fix, 1e6)
          20, ...       %       * ep(3) = ground station G/T, dB                                  (ALICINO fix, 30)
         600, ...       %       * ep(4) = range, km                                               (ALICINO fix, 640)
           5, ...       %       * ep(5) = elevation angle, deg                                    (ALICINO fix, 10)
           5, ...       %       * ep(6) = pointing accuracy, deg                                  (ALICINO fix, 5)
          10, ...       %       * ep(7) = ground station antenna diameter, m
           5, ...       %       * ep(8) = ground station access time, min                         (ALICINO fix, 10)
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