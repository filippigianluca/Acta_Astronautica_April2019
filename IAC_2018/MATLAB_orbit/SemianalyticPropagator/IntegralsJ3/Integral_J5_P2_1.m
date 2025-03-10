function I = Integral_J5_P2_1(Eq0, L0, L)


P1 = Eq0(2);
P2 = Eq0(3);
Q1 = Eq0(4);
Q2 = Eq0(5);

G = 1 + Q1^2 + Q2^2;


LL = [L0 L];


I_temp = LL * ( ...
   (-15/4).*G.^(-4).*(G.^4.*((-21).*P1.^3.*P2.*Q1+(-21).*P1.*P2.*(2+ ...
  P2.^2).*Q1+18.*P1.^4.*Q2+P1.^2.*(41+15.*P2.^2).*Q2+(-1).*((-4)+ ...
  P2.^2+3.*P2.^4).*Q2)+(-28).*G.^2.*((-7).*P1.^3.*P2.*Q1.*(Q1.^2+3.* ...
  Q2.^2)+(-1).*((-1)+P2.^2).*Q2.*((2+3.*P2.^2).*Q1.^2+(2+P2.^2).* ...
  Q2.^2)+P1.^2.*Q2.*((19+21.*P2.^2).*Q1.^2+3.*(7+P2.^2).*Q2.^2)+(-1) ...
  .*P1.*P2.*Q1.*((20+13.*P2.^2).*Q1.^2+3.*(8+P2.^2).*Q2.^2)+2.* ...
  P1.^4.*(3.*Q1.^2.*Q2+5.*Q2.^3))+21.*(Q1.^2+Q2.^2).*((-21).*P1.^3.* ...
  P2.*Q1.*(Q1.^2+5.*Q2.^2)+(-3).*P1.*P2.*Q1.*((26+19.*P2.^2).*Q1.^2+ ...
  (-1).*((-34)+P2.^2).*Q2.^2)+(-1).*((-1)+P2.^2).*Q2.*((8+15.*P2.^2) ...
  .*Q1.^2+(8+3.*P2.^2).*Q2.^2)+P1.^2.*Q2.*((73+111.*P2.^2).*Q1.^2+( ...
  85+3.*P2.^2).*Q2.^2)+6.*P1.^4.*(3.*Q1.^2.*Q2+7.*Q2.^3)))...
    ) + ...
cos(LL) * ( ...
(3/128).*G.^(-4).*((-5).*G.^4.*(43.*(16+48.*P1.^2+5.*P1.^4).*P2.* ...
  Q1+86.*(8+3.*P1.^2).*P2.^3.*Q1+43.*P2.^5.*Q1+(-1).*P1.*(1296+ ...
  2288.*P1.^2+205.*P1.^4).*Q2+(-2).*P1.*(456+119.*P1.^2).*P2.^2.*Q2+ ...
  (-33).*P1.*P2.^4.*Q2)+140.*G.^2.*((-2).*P1.*P2.^2.*Q2.*((678+133.* ...
  P1.^2).*Q1.^2+(78+35.*P1.^2).*Q2.^2)+2.*P2.^3.*Q1.*((218+63.* ...
  P1.^2).*Q1.^2+(34+69.*P1.^2).*Q2.^2)+(-1).*P1.*Q2.*((624+796.* ...
  P1.^2+53.*P1.^4).*Q1.^2+(656+1260.*P1.^2+119.*P1.^4).*Q2.^2)+P2.* ...
  Q1.*((336+708.*P1.^2+55.*P1.^4).*Q1.^2+(368+2004.*P1.^2+265.* ...
  P1.^4).*Q2.^2)+P2.^5.*(31.*Q1.^3+(-7).*Q1.*Q2.^2)+P1.*P2.^4.*(( ...
  -93).*Q1.^2.*Q2+9.*Q2.^3))+42.*(5.*P1.*P2.^4.*Q2.*(325.*Q1.^4+94.* ...
  Q1.^2.*Q2.^2+(-31).*Q2.^4)+(-10).*P2.^3.*Q1.*((484+113.*P1.^2).* ...
  Q1.^4+2.*(196+191.*P1.^2).*Q1.^2.*Q2.^2+5.*((-12)+17.*P1.^2).* ...
  Q2.^4)+10.*P1.*P2.^2.*Q2.*((1836+287.*P1.^2).*Q1.^4+2.*(876+245.* ...
  P1.^2).*Q1.^2.*Q2.^2+(12+35.*P1.^2).*Q2.^4)+(-5).*P2.*Q1.*((664+ ...
  1080.*P1.^2+67.*P1.^4).*Q1.^4+2.*(712+3096.*P1.^2+325.*P1.^4).* ...
  Q1.^2.*Q2.^2+5.*(152+984.*P1.^2+147.*P1.^4).*Q2.^4)+P1.*Q2.*(5.*( ...
  1224+1224.*P1.^2+65.*P1.^4).*Q1.^4+10.*(1272+1960.*P1.^2+147.* ...
  P1.^4).*Q1.^2.*Q2.^2+(6600+13160.*P1.^2+1281.*P1.^4).*Q2.^4)+ ...
  P2.^5.*((-371).*Q1.^5+(-10).*Q1.^3.*Q2.^2+145.*Q1.*Q2.^4)))...
) + ...
cos(2*LL) * ( ...
(-15/32).*G.^(-4).*(G.^4.*((32+45.*P1.^4+168.*P2.^2+53.*P2.^4+2.* ...
  P1.^2.*(76+93.*P2.^2)).*Q1+(-4).*P1.*P2.*(40+31.*P1.^2+9.*P2.^2).* ...
  Q2)+(-56).*G.^2.*((8+6.*P1.^4+53.*P2.^2+19.*P2.^4+9.*P1.^2.*(3+5.* ...
  P2.^2)).*Q1.^3+(-1).*P1.*P2.*(106+61.*P1.^2+45.*P2.^2).*Q1.^2.*Q2+ ...
  (8+27.*P1.^4+9.*P2.^2+(-4).*P2.^4+P1.^2.*(71+51.*P2.^2)).*Q1.* ...
  Q2.^2+3.*P1.*P2.*((-6)+(-7).*P1.^2+P2.^2).*Q2.^3)+42.*((32+15.* ...
  P1.^4+236.*P2.^2+91.*P2.^4+6.*P1.^2.*(14+27.*P2.^2)).*Q1.^5+(-24) ...
  .*P1.*P2.*(24+11.*P1.^2+13.*P2.^2).*Q1.^4.*Q2+2.*(32+69.*P1.^4+ ...
  92.*P2.^2+P2.^4+6.*P1.^2.*(38+45.*P2.^2)).*Q1.^3.*Q2.^2+(-32).* ...
  P1.*P2.*(17+14.*P1.^2+3.*P2.^2).*Q1.^2.*Q2.^3+(32+147.*P1.^4+(-20) ...
  .*P2.^2+(-33).*P2.^4+2.*P1.^2.*(170+69.*P2.^2)).*Q1.*Q2.^4+(-8).* ...
  P1.*P2.*(4+7.*P1.^2+(-3).*P2.^2).*Q2.^5))...
) + ...
cos(3*LL) * ( ...
(-5/128).*G.^(-4).*(G.^4.*(3.*(208+264.*P1.^2+5.*P1.^4).*P2.*Q1+ ...
  6.*(164+39.*P1.^2).*P2.^3.*Q1+75.*P2.^5.*Q1+P1.*(624+1256.*P1.^2+ ...
  127.*P1.^4).*Q2+(-2).*P1.*(12+7.*P1.^2).*P2.^2.*Q2+3.*P1.*P2.^4.* ...
  Q2)+(-56).*G.^2.*(2.*P1.*P2.^2.*Q2.*((-1).*(183+28.*P1.^2).*Q1.^2+ ...
  (57+7.*P1.^2).*Q2.^2)+P2.*Q1.*(2.*(122+165.*P1.^2+8.*P1.^4).* ...
  Q1.^2+(-3).*(36+66.*P1.^2+11.*P1.^4).*Q2.^2)+P1.*Q2.*(((-84)+(-46) ...
  .*P1.^2+P1.^4).*Q1.^2+2.*(118+217.*P1.^2+21.*P1.^4).*Q2.^2)+6.* ...
  P2.^3.*((63+13.*P1.^2).*Q1.^3+(-25).*Q1.*Q2.^2)+15.*P2.^5.*(2.* ...
  Q1.^3+(-1).*Q1.*Q2.^2)+P1.*P2.^4.*((-33).*Q1.^2.*Q2+12.*Q2.^3))+( ...
  -42).*((-3).*P1.*P2.^4.*Q2.*((-119).*Q1.^4+62.*Q1.^2.*Q2.^2+13.* ...
  Q2.^4)+(-2).*P2.^3.*Q1.*((944+169.*P1.^2).*Q1.^4+2.*((-184)+91.* ...
  P1.^2).*Q1.^2.*Q2.^2+(-13).*(32+7.*P1.^2).*Q2.^4)+(-2).*P1.* ...
  P2.^2.*Q2.*((-1).*(1752+227.*P1.^2).*Q1.^4+6.*(96+P1.^2).*Q1.^2.* ...
  Q2.^2+3.*(72+11.*P1.^2).*Q2.^4)+(-1).*P1.*Q2.*((-1).*(936+704.* ...
  P1.^2+25.*P1.^4).*Q1.^4+6.*(88+112.*P1.^2+11.*P1.^4).*Q1.^2.* ...
  Q2.^2+3.*(360+672.*P1.^2+65.*P1.^4).*Q2.^4)+(-3).*P2.*Q1.*((392+ ...
  464.*P1.^2+21.*P1.^4).*Q1.^4+2.*((-8)+320.*P1.^2+23.*P1.^4).* ...
  Q1.^2.*Q2.^2+(-1).*(280+848.*P1.^2+111.*P1.^4).*Q2.^4)+P2.^5.*(( ...
  -155).*Q1.^5+110.*Q1.^3.*Q2.^2+65.*Q1.*Q2.^4)))...
) + ...
cos(4*LL) * ( ...
 (3/16).*G.^(-4).*(5.*G.^4.*(9.*P1.^4.*Q1+P1.^2.*(19+3.*P2.^2).*Q1+ ...
  (-1).*P2.^2.*(19+10.*P2.^2).*Q1+(-29).*P1.^3.*P2.*Q2+(-1).*P1.* ...
  P2.*(38+9.*P2.^2).*Q2)+(-140).*G.^2.*(((-2)+P1.^4+(-20).*P2.^2+( ...
  -6).*P1.^2.*P2.^2+(-9).*P2.^4).*Q1.^3+(-2).*P1.*P2.*(11+8.*P1.^2+ ...
  3.*P2.^2).*Q1.^2.*Q2+(6+15.*P1.^4+22.*P2.^2+7.*P2.^4+P1.^2.*(38+ ...
  24.*P2.^2)).*Q1.*Q2.^2+(-2).*P1.*P2.*(9+7.*P1.^2+2.*P2.^2).*Q2.^3) ...
  +42.*((-1).*(32+285.*P2.^2+125.*P2.^4+35.*P1.^2.*(1+3.*P2.^2)).* ...
  Q1.^5+(-5).*P1.*P2.*(2+11.*P1.^2+(-9).*P2.^2).*Q1.^4.*Q2+10.*(8+ ...
  12.*P1.^4+45.*P2.^2+17.*P2.^4+P1.^2.*(35+33.*P2.^2)).*Q1.^3.* ...
  Q2.^2+(-10).*P1.*P2.*(86+53.*P1.^2+33.*P2.^2).*Q1.^2.*Q2.^3+5.*( ...
  16+48.*P1.^4+43.*P2.^2+11.*P2.^4+9.*P1.^2.*(13+7.*P2.^2)).*Q1.* ...
  Q2.^4+(-5).*P1.*P2.*(26+23.*P1.^2+3.*P2.^2).*Q2.^5))...
) + ...
cos(5 * LL) * ( ...
  (3/128).*G.^(-4).*(G.^4.*(101.*P1.^4.*P2.*Q1+(-1).*P2.^3.*(296+ ...
  39.*P2.^2).*Q1+P1.^2.*(888.*P2.*Q1+94.*P2.^3.*Q1)+45.*P1.^5.*Q2+ ...
  2.*P1.^3.*(148+(-77).*P2.^2).*Q2+(-1).*P1.*P2.^2.*(888+71.*P2.^2) ...
  .*Q2)+56.*G.^2.*(P1.*P2.^4.*Q2.*(53.*Q1.^2+6.*Q2.^2)+2.*P1.* ...
  P2.^2.*Q2.*((309+50.*P1.^2).*Q1.^2+9.*(5+P1.^2).*Q2.^2)+(-2).* ...
  P2.^3.*Q1.*(((-97)+3.*P1.^2).*Q1.^2+(143+38.*P1.^2).*Q2.^2)+P1.* ...
  Q2.*((228+250.*P1.^2+15.*P1.^4).*Q1.^2+(-2).*(38+91.*P1.^2+10.* ...
  P1.^4).*Q2.^2)+(-1).*P2.*Q1.*(2.*((-38)+63.*P1.^2+9.*P1.^4).* ...
  Q1.^2+(228+510.*P1.^2+47.*P1.^4).*Q2.^2)+P2.^5.*(20.*Q1.^3+(-21).* ...
  Q1.*Q2.^2))+(-21).*(P1.*P2.^4.*Q2.*(445.*Q1.^4+806.*Q1.^2.*Q2.^2+( ...
  -23).*Q2.^4)+2.*P1.*P2.^2.*Q2.*((2544+403.*P1.^2).*Q1.^4+2.*(2400+ ...
  397.*P1.^2).*Q1.^2.*Q2.^2+((-48)+7.*P1.^2).*Q2.^4)+2.*P2.^3.*Q1.*( ...
  5.*(256+11.*P1.^2).*Q1.^4+(-2).*(1744+419.*P1.^2).*Q1.^2.*Q2.^2+( ...
  -1).*(544+189.*P1.^2).*Q2.^4)+P1.*Q2.*((1776+1856.*P1.^2+105.* ...
  P1.^4).*Q1.^4+2.*(1872+2144.*P1.^2+135.*P1.^4).*Q1.^2.*Q2.^2+(-1) ...
  .*(1104+2176.*P1.^2+219.*P1.^4).*Q2.^4)+(-1).*P2.*Q1.*(((-1264)+ ...
  96.*P1.^2+71.*P1.^4).*Q1.^4+2.*(2672+5568.*P1.^2+509.*P1.^4).* ...
  Q1.^2.*Q2.^2+(976+2592.*P1.^2+243.*P1.^4).*Q2.^4)+P2.^5.*(245.* ...
  Q1.^5+(-530).*Q1.^3.*Q2.^2+(-71).*Q1.*Q2.^4)))...
) + ...
cos(6 * LL) * ( ...
 (1/32).*G.^(-4).*((-45).*G.^4.*(P1.^4.*Q1+(-6).*P1.^2.*P2.^2.*Q1+ ...
  P2.^4.*Q1+(-4).*P1.^3.*P2.*Q2+4.*P1.*P2.^3.*Q2)+21.*(((-64)+75.* ...
  P1.^4+(-920).*P2.^2+(-525).*P2.^4+10.*P1.^2.*(28+39.*P2.^2)).* ...
  Q1.^5+(-60).*P1.*P2.*(64+23.*P1.^2+41.*P2.^2).*Q1.^4.*Q2+10.*(64+ ...
  21.*P1.^4+488.*P2.^2+189.*P2.^4+2.*P1.^2.*(76+165.*P2.^2)).* ...
  Q1.^3.*Q2.^2+120.*P1.*P2.*((-8)+3.*P1.^2+(-11).*P2.^2).*Q1.^2.* ...
  Q2.^3+(-5).*(64+261.*P1.^4+56.*P2.^2+(-3).*P2.^4+2.*P1.^2.*(292+ ...
  93.*P2.^2)).*Q1.*Q2.^4+60.*P1.*P2.*(16+13.*P1.^2+3.*P2.^2).*Q2.^5) ...
  +(-280).*G.^2.*(3.*P1.^2.*Q1.*((3+5.*P2.^2).*Q1.^2+3.*((-3)+P2.^2) ...
  .*Q2.^2)+P1.*P2.*Q2.*((-3).*(18+13.*P2.^2).*Q1.^2+(18+P2.^2).* ...
  Q2.^2)+P2.^2.*Q1.*((-1).*(9+7.*P2.^2).*Q1.^2+3.*(9+4.*P2.^2).* ...
  Q2.^2)+P1.^4.*(2.*Q1.^3+(-15).*Q1.*Q2.^2)+P1.^3.*P2.*((-15).* ...
  Q1.^2.*Q2+17.*Q2.^3)))...
) + ...
cos(7 * LL) * ( ...
 (-15/128).*G.^(-4).*(G.^4.*(5.*P1.^4.*P2.*Q1+(-10).*P1.^2.*P2.^3.* ...
  Q1+P2.^5.*Q1+P1.^5.*Q2+(-10).*P1.^3.*P2.^2.*Q2+5.*P1.*P2.^4.*Q2)+ ...
  2.*G.^2.*(P1.^4.*P2.*Q1.*(11.*Q1.^2+(-173).*Q2.^2)+(-2).*P1.^3.* ...
  Q2.*(((-204)+P2.^2).*Q1.^2+(68+(-47).*P2.^2).*Q2.^2)+P1.*P2.^2.* ...
  Q2.*((-1).*(1224+203.*P2.^2).*Q1.^2+3.*(136+7.*P2.^2).*Q2.^2)+2.* ...
  P1.^2.*P2.*Q1.*(3.*(68+19.*P2.^2).*Q1.^2+(-1).*(612+31.*P2.^2).* ...
  Q2.^2)+P2.^3.*Q1.*((-1).*(136+25.*P2.^2).*Q1.^2+(408+47.*P2.^2).* ...
  Q2.^2)+P1.^5.*(41.*Q1.^2.*Q2+(-23).*Q2.^3))+(-3).*(P1.*P2.^4.*Q2.* ...
  ((-829).*Q1.^4+34.*Q1.^2.*Q2.^2+47.*Q2.^4)+2.*P1.*P2.^2.*Q2.*((-1) ...
  .*(3156+223.*P1.^2).*Q1.^4+6.*(236+73.*P1.^2).*Q1.^2.*Q2.^2+3.*( ...
  116+23.*P1.^2).*Q2.^4)+P1.*Q2.*(((-1200)+(-296).*P1.^2+15.*P1.^4) ...
  .*Q1.^4+2.*(1200+1928.*P1.^2+149.*P1.^4).*Q1.^2.*Q2.^2+(-1).*(240+ ...
  712.*P1.^2+85.*P1.^4).*Q2.^4)+P2.*Q1.*(((-240)+1032.*P1.^2+59.* ...
  P1.^4).*Q1.^4+(-2).*((-1200)+264.*P1.^2+163.*P1.^4).*Q1.^2.*Q2.^2+ ...
  (-1).*(1200+4632.*P1.^2+529.*P1.^4).*Q2.^4)+P2.^5.*((-105).*Q1.^5+ ...
  450.*Q1.^3.*Q2.^2+(-37).*Q1.*Q2.^4)+P2.^3.*(((-824)+226.*P1.^2).* ...
  Q1.^5+4.*(1244+119.*P1.^2).*Q1.^3.*Q2.^2+(-2).*(428+243.*P1.^2).* ...
  Q1.*Q2.^4)))...
) + ...
cos(8 * LL) * ( ...
(105/32).*G.^(-4).*(2.*G.^2.*((-6).*P1.^2.*P2.^2.*Q1.*(Q1.^2+(-3) ...
  .*Q2.^2)+P2.^4.*Q1.*(Q1.^2+(-3).*Q2.^2)+4.*P1.^3.*P2.*Q2.*((-3).* ...
  Q1.^2+Q2.^2)+(-4).*P1.*P2.^3.*Q2.*((-3).*Q1.^2+Q2.^2)+P1.^4.*( ...
  Q1.^3+(-3).*Q1.*Q2.^2))+(-3).*(24.*P1.^4.*Q1.*Q2.^2.*(Q1.^2+(-1).* ...
  Q2.^2)+P1.^3.*P2.*Q2.*((-11).*Q1.^4+(-74).*Q1.^2.*Q2.^2+17.*Q2.^4) ...
  +P1.*P2.*Q2.*((70+81.*P2.^2).*Q1.^4+(-2).*(70+33.*P2.^2).*Q1.^2.* ...
  Q2.^2+(14+(-3).*P2.^2).*Q2.^4)+P2.^2.*Q1.*(7.*(1+P2.^2).*Q1.^4+( ...
  -2).*(35+23.*P2.^2).*Q1.^2.*Q2.^2+(35+11.*P2.^2).*Q2.^4)+P1.^2.*(( ...
  -7).*(1+3.*P2.^2).*Q1.^5+2.*(35+33.*P2.^2).*Q1.^3.*Q2.^2+((-35)+ ...
  39.*P2.^2).*Q1.*Q2.^4)))...
) + ...
cos(9 * LL) * ( ...
(7/128).*G.^(-4).*((-75).*P1.^4.*P2.*Q1.*(Q1.^4+14.*Q1.^2.*Q2.^2+( ...
  -19).*Q2.^4)+10.*P1.^3.*Q2.*(5.*(52+23.*P2.^2).*Q1.^4+130.*((-4)+ ...
  P2.^2).*Q1.^2.*Q2.^2+(52+(-49).*P2.^2).*Q2.^4)+10.*P1.^2.*P2.*Q1.* ...
  ((156+67.*P2.^2).*Q1.^4+(-10).*(156+31.*P2.^2).*Q1.^2.*Q2.^2+5.*( ...
  156+(-5).*P2.^2).*Q2.^4)+(-15).*P1.*P2.^2.*Q2.*(5.*(104+25.*P2.^2) ...
  .*Q1.^4+(-130).*(8+P2.^2).*Q1.^2.*Q2.^2+(104+P2.^2).*Q2.^4)+(-1).* ...
  P2.^3.*Q1.*((520+119.*P2.^2).*Q1.^4+(-10).*(520+83.*P2.^2).* ...
  Q1.^2.*Q2.^2+5.*(520+47.*P2.^2).*Q2.^4)+P1.^5.*(145.*Q1.^4.*Q2+( ...
  -650).*Q1.^2.*Q2.^3+101.*Q2.^5)+10.*G.^2.*(5.*P1.^4.*P2.*Q1.*( ...
  Q1.^2+(-3).*Q2.^2)+(-10).*P1.^2.*P2.^3.*Q1.*(Q1.^2+(-3).*Q2.^2)+ ...
  P2.^5.*Q1.*(Q1.^2+(-3).*Q2.^2)+10.*P1.^3.*P2.^2.*Q2.*((-3).*Q1.^2+ ...
  Q2.^2)+(-5).*P1.*P2.^4.*Q2.*((-3).*Q1.^2+Q2.^2)+P1.^5.*(3.*Q1.^2.* ...
  Q2+(-1).*Q2.^3)))...
) + ...
cos(10 * LL) * ( ...
(-189/32).*G.^(-4).*((-4).*P1.^3.*P2.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.* ...
  Q2.^2+Q2.^4)+4.*P1.*P2.^3.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.*Q2.^2+ ...
  Q2.^4)+(-6).*P1.^2.*P2.^2.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.* ...
  Q2.^4)+P2.^4.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+P1.^4.*( ...
  Q1.^5+(-10).*Q1.^3.*Q2.^2+5.*Q1.*Q2.^4))...
) + ...
cos(11 * LL) * ( ...
(-63/128).*G.^(-4).*((-10).*P1.^3.*P2.^2.*Q2.*(5.*Q1.^4+(-10).* ...
  Q1.^2.*Q2.^2+Q2.^4)+5.*P1.*P2.^4.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.* ...
  Q2.^2+Q2.^4)+5.*P1.^4.*P2.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.* ...
  Q2.^4)+(-10).*P1.^2.*P2.^3.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.* ...
  Q2.^4)+P2.^5.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+P1.^5.*(5.* ...
  Q1.^4.*Q2+(-10).*Q1.^2.*Q2.^3+Q2.^5))...
) + ...
 sin(LL) * ( ...
(15/128).*G.^(-4).*(G.^4.*(35.*P1.^5.*Q1+2.*P1.^3.*(296+121.* ...
  P2.^2).*Q1+3.*P1.*(208+656.*P2.^2+69.*P2.^4).*Q1+(-127).*P1.^4.* ...
  P2.*Q2+(-2).*P1.^2.*P2.*(504+41.*P2.^2).*Q2+P2.*(16+368.*P2.^2+ ...
  45.*P2.^4).*Q2)+(-28).*G.^2.*((-1).*P1.^4.*P2.*Q2.*(119.*Q1.^2+ ...
  45.*Q2.^2)+(-2).*P1.^2.*P2.*Q2.*(7.*(90+13.*P2.^2).*Q1.^2+(-3).*(( ...
  -42)+P2.^2).*Q2.^2)+2.*P1.^3.*Q1.*(3.*(34+19.*P2.^2).*Q1.^2+(286+ ...
  71.*P2.^2).*Q2.^2)+P1.*Q1.*((304+1212.*P2.^2+145.*P2.^4).*Q1.^2+ ...
  3.*(112+100.*P2.^2+(-7).*P2.^4).*Q2.^2)+P2.*Q2.*(((-16)+388.* ...
  P2.^2+57.*P2.^4).*Q1.^2+(16+116.*P2.^2+11.*P2.^4).*Q2.^2)+P1.^5.*( ...
  9.*Q1.^3+43.*Q1.*Q2.^2))+42.*((-1).*P1.^4.*P2.*Q2.*(251.*Q1.^4+ ...
  450.*Q1.^2.*Q2.^2+63.*Q2.^4)+(-2).*P1.^2.*P2.*Q2.*((1668+305.* ...
  P2.^2).*Q1.^4+2.*(852+59.*P2.^2).*Q1.^2.*Q2.^2+(132+(-19).*P2.^2) ...
  .*Q2.^4)+2.*P1.^3.*Q1.*((156+101.*P2.^2).*Q1.^4+2.*(444+179.* ...
  P2.^2).*Q1.^2.*Q2.^2+35.*(20+3.*P2.^2).*Q2.^4)+P1.*Q1.*((600+ ...
  2664.*P2.^2+343.*P2.^4).*Q1.^4+2.*(648+1224.*P2.^2+25.*P2.^4).* ...
  Q1.^2.*Q2.^2+(696+(-24).*P2.^2+(-109).*P2.^4).*Q2.^4)+P2.*Q2.*((( ...
  -56)+1000.*P2.^2+161.*P2.^4).*Q1.^4+2.*((-8)+552.*P2.^2+67.*P2.^4) ...
  .*Q1.^2.*Q2.^2+(40+168.*P2.^2+13.*P2.^4).*Q2.^4)+P1.^5.*(11.* ...
  Q1.^5+106.*Q1.^3.*Q2.^2+119.*Q1.*Q2.^4))) ...
) + ...
sin(2*LL) * ( ...
(15/32).*G.^(-4).*(G.^4.*((-52).*P1.^3.*P2.*Q1+4.*P1.*P2.*((-4)+ ...
  9.*P2.^2).*Q1+111.*P1.^4.*Q2+6.*P1.^2.*(40+9.*P2.^2).*Q2+(32+80.* ...
  P2.^2+31.*P2.^4).*Q2)+(-28).*G.^2.*((-2).*P1.^3.*P2.*Q1.*(Q1.^2+ ...
  49.*Q2.^2)+2.*P1.*P2.*Q1.*((18+19.*P2.^2).*Q1.^2+(-7).*(10+3.* ...
  P2.^2).*Q2.^2)+2.*P1.^2.*Q2.*(9.*(3+2.*P2.^2).*Q1.^2+(71+12.* ...
  P2.^2).*Q2.^2)+Q2.*((16+106.*P2.^2+47.*P2.^4).*Q1.^2+(16+18.* ...
  P2.^2+5.*P2.^4).*Q2.^2)+P1.^4.*(21.*Q1.^2.*Q2+67.*Q2.^3))+42.*(4.* ...
  P1.^3.*P2.*Q1.*(3.*Q1.^4+(-42).*Q1.^2.*Q2.^2+(-77).*Q2.^4)+4.*P1.* ...
  P2.*Q1.*(3.*(10+9.*P2.^2).*Q1.^4+(-42).*(2+P2.^2).*Q1.^2.*Q2.^2+( ...
  -7).*(14+3.*P2.^2).*Q2.^4)+Q2.*((32+288.*P2.^2+143.*P2.^4).*Q1.^4+ ...
  2.*(32+136.*P2.^2+45.*P2.^4).*Q1.^2.*Q2.^2+(32+16.*P2.^2+3.*P2.^4) ...
  .*Q2.^4)+3.*P1.^4.*(5.*Q1.^4.*Q2+46.*Q1.^2.*Q2.^3+49.*Q2.^5)+ ...
  P1.^2.*(2.*(16+3.*P2.^2).*Q1.^4.*Q2+92.*(4+3.*P2.^2).*Q1.^2.* ...
  Q2.^3+2.*(152+15.*P2.^2).*Q2.^5))) ...
) + ...
sin(3 * LL) * ( ...
(-5/128).*G.^(-4).*(G.^4.*(63.*P1.^5.*Q1+6.*P1.^3.*(148+43.*P2.^2) ...
  .*Q1+3.*P1.*(208+360.*P2.^2+17.*P2.^4).*Q1+(-243).*P1.^4.*P2.*Q2+( ...
  -2).*P1.^2.*P2.*(948+73.*P2.^2).*Q2+(-1).*P2.*(624+616.*P2.^2+47.* ...
  P2.^4).*Q2)+(-56).*G.^2.*((-3).*P1.^4.*P2.*Q2.*(30.*Q1.^2+17.* ...
  Q2.^2)+(-2).*P1.^2.*P2.*Q2.*(5.*(87+11.*P2.^2).*Q1.^2+3.*(57+2.* ...
  P2.^2).*Q2.^2)+6.*P1.^3.*Q1.*((17+7.*P2.^2).*Q1.^2+(97+22.*P2.^2) ...
  .*Q2.^2)+(-1).*P2.*Q2.*((420+550.*P2.^2+44.*P2.^4).*Q1.^2+(68+22.* ...
  P2.^2+P2.^4).*Q2.^2)+P1.*Q1.*(2.*(38+75.*P2.^2+2.*P2.^4).*Q1.^2+ ...
  3.*(132+210.*P2.^2+13.*P2.^4).*Q2.^2)+P1.^5.*(6.*Q1.^3+45.*Q1.* ...
  Q2.^2))+42.*((-3).*P1.^4.*P2.*Q2.*(97.*Q1.^4+286.*Q1.^2.*Q2.^2+ ...
  53.*Q2.^4)+(-2).*P1.^2.*P2.*Q2.*(9.*(184+29.*P2.^2).*Q1.^4+2.*( ...
  1824+179.*P2.^2).*Q1.^2.*Q2.^2+(456+(-7).*P2.^2).*Q2.^4)+2.* ...
  P1.^3.*Q1.*((104+49.*P2.^2).*Q1.^4+2.*(704+259.*P2.^2).*Q1.^2.* ...
  Q2.^2+(1624+269.*P2.^2).*Q2.^4)+P2.*Q2.*((-3).*(776+1184.*P2.^2+ ...
  101.*P2.^4).*Q1.^4+(-2).*(1032+848.*P2.^2+49.*P2.^4).*Q1.^2.* ...
  Q2.^2+((-120)+64.*P2.^2+5.*P2.^4).*Q2.^4)+3.*P1.*Q1.*((40+32.* ...
  P2.^2+(-11).*P2.^4).*Q1.^4+2.*(408+1040.*P2.^2+87.*P2.^4).*Q1.^2.* ...
  Q2.^2+(648+640.*P2.^2+17.*P2.^4).*Q2.^4)+P1.^5.*(11.*Q1.^5+178.* ...
  Q1.^3.*Q2.^2+271.*Q1.*Q2.^4)))...
) + ...
sin(4 * LL) * ( ...
  (-3/16).*G.^(-4).*(5.*G.^4.*(17.*P1.^3.*P2.*Q1+P1.*P2.*(38+21.* ...
  P2.^2).*Q1+12.*P1.^4.*Q2+P1.^2.*(19+(-15).*P2.^2).*Q2+(-1).* ...
  P2.^2.*(19+7.*P2.^2).*Q2)+70.*G.^2.*((-2).*P1.^3.*P2.*Q1.*(9.* ...
  Q1.^2+7.*Q2.^2)+(-2).*P1.*P2.*Q1.*((20+11.*P2.^2).*Q1.^2+(16+9.* ...
  P2.^2).*Q2.^2)+Q2.*((12+82.*P2.^2+31.*P2.^4).*Q1.^2+(-1).*(4+2.* ...
  P2.^2+P2.^4).*Q2.^2)+P1.^2.*(2.*(19+30.*P2.^2).*Q1.^2.*Q2+(-38).* ...
  Q2.^3)+P1.^4.*(9.*Q1.^2.*Q2+(-19).*Q2.^3))+(-84).*((-5).*P1.^3.* ...
  P2.*Q1.*(7.*Q1.^4+38.*Q1.^2.*Q2.^2+(-5).*Q2.^4)+30.*P1.^4.*Q2.*( ...
  Q1.^4+Q1.^2.*Q2.^2+(-2).*Q2.^4)+(-5).*P1.*P2.*Q1.*((16+9.*P2.^2).* ...
  Q1.^4+2.*(40+21.*P2.^2).*Q1.^2.*Q2.^2+(-1).*(8+3.*P2.^2).*Q2.^4)+ ...
  5.*P1.^2.*Q2.*(9.*(3+5.*P2.^2).*Q1.^4+2.*(11+15.*P2.^2).*Q1.^2.* ...
  Q2.^2+(-1).*(25+3.*P2.^2).*Q2.^4)+Q2.*(5.*(10+73.*P2.^2+29.*P2.^4) ...
  .*Q1.^4+10.*(2+9.*P2.^2+2.*P2.^4).*Q1.^2.*Q2.^2+(-1).*(14+15.* ...
  P2.^2+5.*P2.^4).*Q2.^4)))...
) + ...
sin(5 * LL) * ( ...
 (3/128).*G.^(-4).*(G.^4.*(35.*P1.^5.*Q1+P1.^3.*(296+(-54).*P2.^2) ...
  .*Q1+(-1).*P1.*P2.^2.*(888+121.*P2.^2).*Q1+(-151).*P1.^4.*P2.*Q2+ ...
  6.*P1.^2.*P2.*((-148)+P2.^2).*Q2+P2.^3.*(296+29.*P2.^2).*Q2)+(-56) ...
  .*G.^2.*(35.*P1.^5.*Q1.*Q2.^2+(-1).*P1.^4.*P2.*Q2.*(16.*Q1.^2+45.* ...
  Q2.^2)+2.*P1.^2.*P2.*Q2.*(3.*(11+9.*P2.^2).*Q1.^2+(-1).*(159+8.* ...
  P2.^2).*Q2.^2)+P2.*Q2.*((228+434.*P2.^2+38.*P2.^4).*Q1.^2+(-1).*( ...
  76+46.*P2.^2+3.*P2.^4).*Q2.^2)+P1.*Q1.*((-2).*(38+177.*P2.^2+21.* ...
  P2.^4).*Q1.^2+(228+174.*P2.^2+5.*P2.^4).*Q2.^2)+P1.^3.*((-34).*(1+ ...
  P2.^2).*Q1.^3+2.*(199+24.*P2.^2).*Q1.*Q2.^2))+(-21).*(P1.^4.*P2.* ...
  Q2.*((-179).*Q1.^4+870.*Q1.^2.*Q2.^2+345.*Q2.^4)+2.*P1.^2.*P2.* ...
  Q2.*((-1).*(2160+541.*P2.^2).*Q1.^4+2.*(1632+109.*P2.^2).*Q1.^2.* ...
  Q2.^2+5.*(240+11.*P2.^2).*Q2.^4)+P2.*Q2.*((-1).*(4496+7552.*P2.^2+ ...
  647.*P2.^4).*Q1.^4+2.*(848+608.*P2.^2+39.*P2.^4).*Q1.^2.*Q2.^2+( ...
  560+320.*P2.^2+21.*P2.^4).*Q2.^4)+P1.*Q1.*((720+3072.*P2.^2+359.* ...
  P2.^4).*Q1.^4+2.*(48+1632.*P2.^2+221.*P2.^4).*Q1.^2.*Q2.^2+(-1).*( ...
  3696+4416.*P2.^2+301.*P2.^4).*Q2.^4)+P1.^5.*(11.*Q1.^5+(-110).* ...
  Q1.^3.*Q2.^2+(-505).*Q1.*Q2.^4)+P1.^3.*((416+306.*P2.^2).*Q1.^5+ ...
  4.*((-224)+51.*P2.^2).*Q1.^3.*Q2.^2+(-10).*(592+87.*P2.^2).*Q1.* ...
  Q2.^4)))... 
) + ...
sin(6 * LL) * ( ...
 (1/32).*G.^(-4).*(45.*G.^4.*(4.*P1.^3.*P2.*Q1+(-4).*P1.*P2.^3.*Q1+ ...
  P1.^4.*Q2+(-6).*P1.^2.*P2.^2.*Q2+P2.^4.*Q2)+140.*G.^2.*((-2).* ...
  P1.^3.*P2.*Q1.*(Q1.^2+33.*Q2.^2)+6.*P1.^2.*Q2.*((9+6.*P2.^2).* ...
  Q1.^2+((-3)+4.*P2.^2).*Q2.^2)+P2.^2.*Q2.*((-3).*(18+11.*P2.^2).* ...
  Q1.^2+(18+5.*P2.^2).*Q2.^2)+2.*P1.*P2.*Q1.*((18+19.*P2.^2).*Q1.^2+ ...
  (-3).*(18+7.*P2.^2).*Q2.^2)+P1.^4.*(21.*Q1.^2.*Q2+(-13).*Q2.^3))+( ...
  -21).*(60.*P1.^3.*P2.*Q1.*(3.*Q1.^4+(-38).*Q1.^2.*Q2.^2+(-25).* ...
  Q2.^4)+60.*P1.*P2.*Q1.*((20+17.*P2.^2).*Q1.^4+(-2).*(28+9.*P2.^2) ...
  .*Q1.^2.*Q2.^2+(-1).*(44+19.*P2.^2).*Q2.^4)+10.*P1.^2.*Q2.*((32+( ...
  -39).*P2.^2).*Q1.^4+2.*(184+183.*P2.^2).*Q1.^2.*Q2.^2+((-80)+21.* ...
  P2.^2).*Q2.^4)+Q2.*((-5).*(64+704.*P2.^2+339.*P2.^4).*Q1.^4+10.*( ...
  64+272.*P2.^2+75.*P2.^4).*Q1.^2.*Q2.^2+((-64)+160.*P2.^2+45.* ...
  P2.^4).*Q2.^4)+15.*P1.^4.*(15.*Q1.^4.*Q2+82.*Q1.^2.*Q2.^3+(-29).* ...
  Q2.^5)))...
) + ...
sin(7 * LL) * ( ...
 (-15/128).*G.^(-4).*(G.^4.*(P1.^5.*Q1+(-10).*P1.^3.*P2.^2.*Q1+5.* ...
  P1.*P2.^4.*Q1+(-5).*P1.^4.*P2.*Q2+10.*P1.^2.*P2.^3.*Q2+(-1).* ...
  P2.^5.*Q2)+2.*G.^2.*(P2.^3.*Q2.*((408+61.*P2.^2).*Q1.^2+(-1).*( ...
  136+11.*P2.^2).*Q2.^2)+(-2).*P1.^2.*P2.*Q2.*((612+101.*P2.^2).* ...
  Q1.^2+((-204)+13.*P2.^2).*Q2.^2)+2.*P1.^3.*Q1.*((68+23.*P2.^2).* ...
  Q1.^2+((-204)+71.*P2.^2).*Q2.^2)+P1.*P2.^2.*Q1.*((-1).*(408+91.* ...
  P2.^2).*Q1.^2+(1224+133.*P2.^2).*Q2.^2)+P1.^5.*(9.*Q1.^3+(-55).* ...
  Q1.*Q2.^2)+P1.^4.*P2.*((-103).*Q1.^2.*Q2+81.*Q2.^3))+(-3).*( ...
  P1.^4.*P2.*Q2.*((-299).*Q1.^4+(-226).*Q1.^2.*Q2.^2+217.*Q2.^4)+( ...
  -2).*P1.^3.*Q1.*(((-52)+3.*P2.^2).*Q1.^4+(-2).*(556+291.*P2.^2).* ...
  Q1.^2.*Q2.^2+7.*(196+P2.^2).*Q2.^4)+(-2).*P1.^2.*P2.*Q2.*(3.*(452+ ...
  51.*P2.^2).*Q1.^4+2.*(1092+251.*P2.^2).*Q1.^2.*Q2.^2+(-1).*(708+ ...
  19.*P2.^2).*Q2.^4)+P2.*Q2.*((1200+3304.*P2.^2+361.*P2.^4).*Q1.^4+( ...
  -2).*(1200+1672.*P2.^2+117.*P2.^4).*Q1.^2.*Q2.^2+(240+8.*P2.^2+( ...
  -3).*P2.^4).*Q2.^4)+P1.*Q1.*((-1).*(240+1752.*P2.^2+289.*P2.^4).* ...
  Q1.^4+2.*(1200+3864.*P2.^2+353.*P2.^4).*Q1.^2.*Q2.^2+((-1200)+ ...
  1032.*P2.^2+179.*P2.^4).*Q2.^4)+P1.^5.*(11.*Q1.^5+106.*Q1.^3.* ...
  Q2.^2+(-273).*Q1.*Q2.^4)))...
) + ...
sin(8 * LL) * ( ...
(-105/32).*G.^(-4).*(2.*G.^2.*(4.*P1.^3.*P2.*Q1.*(Q1.^2+(-3).* ...
  Q2.^2)+(-4).*P1.*P2.^3.*Q1.*(Q1.^2+(-3).*Q2.^2)+6.*P1.^2.*P2.^2.* ...
  Q2.*((-3).*Q1.^2+Q2.^2)+(-1).*P2.^4.*Q2.*((-3).*Q1.^2+Q2.^2)+ ...
  P1.^4.*(3.*Q1.^2.*Q2+(-1).*Q2.^3))+3.*(6.*P1.^4.*Q2.*(Q1.^4+(-6).* ...
  Q1.^2.*Q2.^2+Q2.^4)+P1.^3.*P2.*Q1.*((-7).*Q1.^4+(-26).*Q1.^2.* ...
  Q2.^2+61.*Q2.^4)+P1.^2.*Q2.*((35+69.*P2.^2).*Q1.^4+2.*((-35)+3.* ...
  P2.^2).*Q1.^2.*Q2.^2+(7+(-15).*P2.^2).*Q2.^4)+(-1).*P2.^2.*Q2.*(( ...
  35+29.*P2.^2).*Q1.^4+(-2).*(35+17.*P2.^2).*Q1.^2.*Q2.^2+(7+P2.^2) ...
  .*Q2.^4)+P1.*P2.*Q1.*(7.*(2+3.*P2.^2).*Q1.^4+(-2).*(70+57.*P2.^2) ...
  .*Q1.^2.*Q2.^2+(70+9.*P2.^2).*Q2.^4))) ...
) + ...
sin(9 * LL) * ( ...
 (7/128).*G.^(-4).*((-75).*P1.^4.*P2.*Q2.*(Q1.^4+(-26).*Q1.^2.* ...
  Q2.^2+5.*Q2.^4)+10.*P1.^3.*Q1.*((52+41.*P2.^2).*Q1.^4+(-10).*(52+ ...
  5.*P2.^2).*Q1.^2.*Q2.^2+5.*(52+(-31).*P2.^2).*Q2.^4)+(-15).*P1.* ...
  P2.^2.*Q1.*((104+31.*P2.^2).*Q1.^4+(-10).*(104+19.*P2.^2).*Q1.^2.* ...
  Q2.^2+5.*(104+7.*P2.^2).*Q2.^4)+10.*P1.^2.*P2.*Q2.*((-5).*(156+ ...
  49.*P2.^2).*Q1.^4+130.*(12+P2.^2).*Q1.^2.*Q2.^2+((-156)+23.*P2.^2) ...
  .*Q2.^4)+P2.^3.*Q2.*(5.*(520+101.*P2.^2).*Q1.^4+(-650).*(8+P2.^2) ...
  .*Q1.^2.*Q2.^2+(520+29.*P2.^2).*Q2.^4)+P1.^5.*(11.*Q1.^5+(-470).* ...
  Q1.^3.*Q2.^2+415.*Q1.*Q2.^4)+10.*G.^2.*((-10).*P1.^3.*P2.^2.*Q1.*( ...
  Q1.^2+(-3).*Q2.^2)+5.*P1.*P2.^4.*Q1.*(Q1.^2+(-3).*Q2.^2)+5.* ...
  P1.^4.*P2.*Q2.*((-3).*Q1.^2+Q2.^2)+(-10).*P1.^2.*P2.^3.*Q2.*((-3) ...
  .*Q1.^2+Q2.^2)+P2.^5.*Q2.*((-3).*Q1.^2+Q2.^2)+P1.^5.*(Q1.^3+(-3).* ...
  Q1.*Q2.^2)))...
) + ...
sin(10 * LL) * ( ...
 (189/32).*G.^(-4).*((-6).*P1.^2.*P2.^2.*Q2.*(5.*Q1.^4+(-10).* ...
  Q1.^2.*Q2.^2+Q2.^4)+P2.^4.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.*Q2.^2+ ...
  Q2.^4)+4.*P1.^3.*P2.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+(-4) ...
  .*P1.*P2.^3.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+P1.^4.*(5.* ...
  Q1.^4.*Q2+(-10).*Q1.^2.*Q2.^3+Q2.^5))...
) + ...
sin(11 * LL) * ( ...
  (-63/128).*G.^(-4).*((-5).*P1.^4.*P2.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.* ...
  Q2.^2+Q2.^4)+10.*P1.^2.*P2.^3.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.*Q2.^2+ ...
  Q2.^4)+(-1).*P2.^5.*Q2.*(5.*Q1.^4+(-10).*Q1.^2.*Q2.^2+Q2.^4)+(-10) ...
  .*P1.^3.*P2.^2.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+5.*P1.* ...
  P2.^4.*Q1.*(Q1.^4+(-10).*Q1.^2.*Q2.^2+5.*Q2.^4)+P1.^5.*(Q1.^5+( ...
  -10).*Q1.^3.*Q2.^2+5.*Q1.*Q2.^4))...
) ...
;




I = I_temp(2:end)  - I_temp(1);