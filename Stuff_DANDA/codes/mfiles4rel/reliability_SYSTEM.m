function rel = reliability_SYSTEM(t, d,u)

%assume variables are just stuck one sequence after the previous one and extract the subsequence for systems
dseq = [14,8]; useq = [13,12];
d_aocs = d(1:14);
u_aocs = u(1:8);
d_power = d(15:27);
u_power = u(9:20);


relAOCS = reliability_AOCS(t, d_aocs, u_aocs);
relPOWER = reliability_POWER(t, d_power, u_power);

rel = relAOCS * relPOWER;

return
