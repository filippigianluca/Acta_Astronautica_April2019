% prova integrale
TM = 1000;



% 1. integrale del prodotto
P1 = integral(@satellite_survival,0,TM);


% 2. prodotto integrali
% integral(@(t) (1-exp(-(t./3831).^0.7182)),0,TM)
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

for i=1:9
P2(i) = integral(@(t) wblcdf(t, scale(i), beta(i)),0,TM);
end

