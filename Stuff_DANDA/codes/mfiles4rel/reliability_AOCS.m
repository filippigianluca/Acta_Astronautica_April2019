function rel = reliability_AOCS(t, d,u)

lambda = 6206945;
k = 0.3375;

validD = [7,9,10];
validU = [2,3,4,10,11,12,13];

noValidD = length(validD);
noValidU = length(validU);

noValid = noValidD + noValidU;

R0 = exp(-(t/lambda)^k);

w = [1.31916614
  0.59614842
 -0.3045549
 -1.75893981
 -0.47402175
  -2.73295655
  0.61283252
  0.02050589
  0.55533248
  1.08875295];


[lbu, ubu, lbd, ubd] = bounds_sensitivity_exact_subsystems_AOCS();

x = zeros(1,noValid);

shif = 0;
for i=1:noValidD
    position = validD(i);
    add = -lbd(position);
    coef = 1./(ubd(position)+add);
    x(i) = coef*(d(position)+add);
    shif = shif + 1;
end

for i=1:noValidU
    position = validU(i);
    add = -lbu(position);
    coef = 1./(ubu(position)+add);
    x(i+shif) = coef*(u(position)+add);
end


exponent = x * w;
rel = real(R0^exp(exponent));

return
