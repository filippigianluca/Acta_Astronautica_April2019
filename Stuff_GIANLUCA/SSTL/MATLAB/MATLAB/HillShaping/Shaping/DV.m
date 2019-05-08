function xdot = DV(t,x,F,r,G,i,F_n)

xdot = abs(F(t)./(G(t)./r(t).^2-r(t).*cot(i(t)).*sin(t).*F_n(t)./G(t))) ;

end