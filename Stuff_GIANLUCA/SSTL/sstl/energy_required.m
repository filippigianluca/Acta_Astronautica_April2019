function [e_req, e_sizing] = energy_required(t_on_off,eff)
eff =1.1;

P_out_on = 1925.0;
P_in_off = (2505.0-1539.0)/eff;
e_req = 0;
e = 0;
e0 = 0;
for i = 1:size(t_on_off,1)
    if i>1
        e0 = max(0.0, e-(t_on_off(i,1)-t_on_off(i-1,2))*P_in_off); % offset if the battery didn't have time to charge completely
        e_in_sizing(i) = (t_on_off(i,1)-t_on_off(i-1,2))*P_in_off; %(Wh)
    end
    e = e0 + (t_on_off(i,2) - t_on_off(i,1))*P_out_on;
    e_req = max(e,e_req);
    
    e_out_sizing(i) = e;
end

e_sizing.in = e_in_sizing;
e_sizing.out = e_out_sizing;
global e0g
e0g = [e0g e0];
return