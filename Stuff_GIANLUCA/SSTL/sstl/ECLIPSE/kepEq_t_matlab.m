function dt = kepEq_t_matlab(th, a, e, mu, th_0, t_0)

%% final position

th = mod(th,2*pi);
th_0 = mod(th_0,2*pi);
% if th>pi && th<2*pi
%     th = th-2*pi
% end

% E = atan2((1-e^2)^.5*sin(th), (e+cos(th)));
% E = acos((e+cos(th))/(1+e*cos(th)));
E = 2*atan2((1-e)^.5*sin(th/2), (1+e)^.5*cos(th/2));

E = mod(E,2*pi);
% if th > 0 && th < pi
%     if E>pi
%         E = E - pi;
%     end
% elseif th > pi && th < 2*pi
%     if E<pi
%         E = E + pi;
%     end
% end
M = E - e*sin(E);


%% initial position

% E_0 = acos((e+cos(th_0))/(1+e*cos(th_0)));
% E_0 = atan2((1-e^2)^.5*sin(th_0), (e+cos(th_0)));
E_0 = 2*atan2((1-e)^.5*sin(th_0/2), (1+e)^.5*cos(th_0/2));

E_0 = mod(E_0,2*pi);
% if th_0 > 0 && th_0 < pi
%     if E_0>pi
%         E_0 = E_0 - pi;
%     end
% elseif th_0 > pi && th_0 < 2*pi
%     if E_0<pi
%         E_0 = E_0 + pi;
%     end
% end

M_0 = E_0 - e*sin(E_0);

n = (mu/a^3)^.5;

t = t_0 + (M-M_0)/n;
% dt = abs((M-M_0)/n);
dt = (M-M_0)/n;

% p = calculateEllipse(0, 0, a, e, 0, th_0, th);


return