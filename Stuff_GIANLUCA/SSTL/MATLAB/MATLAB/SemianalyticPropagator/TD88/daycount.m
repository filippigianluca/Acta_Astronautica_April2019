function d=daycount(t)
% d=datenum(datetime(t))-datenum(datetime(to));
% d=(t-to)/86400;
% Y=year(t);
% year_start=datetime(1,1,Y);
% d=(datenum(t)-datenum(year_start)); % turn to radians

% daycount=(datenum(datetime(t))-datenum(datetime(to)));
% if daycount<=360
%     d=daycount;
% end
% if 360<daycount<=720
%     d=daycount-360;
% end
% 
% if 720<daycount<=1080;
%     d=daycount-2*360;
% end
% 
% if 1080<daycount<1440 
%     d=daycount-3*360;
% end
% 
% if 1440<daycount<1800
%     d=daycount-4*360;
% end
% if 1800<daycount<2160
%     d=daycount-5*360;
% end

d=((t-floor(t))/0.00274); % in radians for sin(d-p3) to make more sense
end
 