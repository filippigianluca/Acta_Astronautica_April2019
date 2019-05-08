function x = acos2(cosx,H)

x = wrapTo2Pi(sign(H).*acos(cosx)) ;

end