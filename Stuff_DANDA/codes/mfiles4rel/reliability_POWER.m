function rel = reliability_POWER(t, d,u)

lambda = 169.272;
k = 0.5021;

validD = [2,5,8];
validU = [7,11];

noValidD = length(validD);
noValidU = length(validU);

noValid = noValidD + noValidU;

R0 = exp(-(t/lambda)^k);

w = [ 0.72494093
   -0.33792166
   -1.61555565
    0.82458707
   -0.02078624];


[lbu, ubu, lbd, ubd] = bounds_sensitivity_exact_subsystems_POWER();

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
