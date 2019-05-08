function out = satellite_survival(t)

beta = [ 0.7182
    0.3375
    1.4560
    0.356
    0.8874
    0.746
    0.5021
    0.4035
    0.3939];

scale = [ 3831
    6206945
    408
    21308746
    7983
    7733
    169.272
    1965868
    400982];
    
survivals = 1. - wblcdf(t, scale, beta);
out = prod(survivals);
end