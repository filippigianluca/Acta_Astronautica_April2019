function I = Integral_J5_P2_2(Eq0, L)


P1 = Eq0(2);
P2 = Eq0(3);
Q1 = Eq0(4);
Q2 = Eq0(5);

G = 1 + Q1^2 + Q2^2;


I = L * ( ...
   ...
    ) + ...
cos(L) * ( ...
...
) + ...
cos(2*L) * ( ...
...
) + ...
cos(3*L) * ( ...
...
) + ...
cos(4*L) * ( ...
 ...
) + ...
cos(5 * L) * ( ...
 ...
) + ...
cos(6 * L) * ( ...
 ...
) + ...
cos(7 * L) * ( ...
...
) + ...
cos(8 * L) * ( ...
...
) + ...
cos(9 * L) * ( ...
...
) + ...
cos(10 * L) * ( ...
...
) + ...
cos(11 * L) * ( ...
...
) + ...
 sin(L) * ( ...
 ...
) + ...
sin(2*L) * ( ...
 ...
) + ...
sin(3 * L) * ( ...
...
) + ...
sin(4 * L) * ( ...
  ...
) + ...
sin(5 * L) * ( ...
 ... 
) + ...
sin(6 * L) * ( ...
 ...
) + ...
sin(7 * L) * ( ...
 ...
) + ...
sin(8 * L) * ( ...
 ...
) + ...
sin(9 * L) * ( ...
 ...
) + ...
sin(10 * L) * ( ...
 ...
) + ...
sin(11 * L) * ( ...
 ...
) ...
;

