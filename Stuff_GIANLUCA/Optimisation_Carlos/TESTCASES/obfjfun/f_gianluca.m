function [f] = f_gianluca(d,u,par)

f1 = 10*u(1)^2 + abs(u(2))*u(5)^2 + u(6)^4/100 + d(1)^2*abs(d(2));
f2 = abs(u(3)) + u(4)^2 * abs(u(5))/10 + u(6)^2 + abs(d(1));
f=f1+f2;

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

global history;
history = [history; d u f];
% figure(1)
% plot(history(:,end),'b.');
% title('f');
% figure(2)
% subplot(1,2,1);
% plot(history(:,1),'k.');
% title('d_1');
% subplot(1,2,2);
% plot(history(:,2),'k.');
% title('d_2');

end
