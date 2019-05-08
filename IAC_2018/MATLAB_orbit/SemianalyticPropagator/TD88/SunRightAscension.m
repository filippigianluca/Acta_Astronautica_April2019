function alpha = SunRightAscension(mjd2000)

T = (mjd2000 ) / 36525;

L0 = mod((280.46645 + 36000.76983 * T +0.0003032 * T^2), 360);

M = mod((357.52910 + 35999.05030 * T - 0.0001559 * T^2 - 0.00000048 * T^3), 360);

C = (1.914600 - 0.004817 * T -0.000014 * T^2) * sind(M) + ...
    (0.019993 - 0.000101 * T) * sind(2*M) + ...
    0.000290 * sind(3*M);

theta = L0 + C;

epsilon = 23.26;

alpha = atan2d(cosd(epsilon) * sind(theta), cosd(theta));

alpha = alpha * pi/180;


end