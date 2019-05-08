function [dvlog] = obfjun_dv_transfer(d,u,par)
global nfevalglobal;
nfevalglobal = nfevalglobal + 1;

	mjd = d(1)+u(1);
	tof = d(2)+u(2);

	ephref0 = par.ephref0 + u(3:8);
	ephref1 = par.ephref1 + u(9:14);
	epoch_ref = par.epoch_ref;

	% origin
	epoch_diff0 = (mjd-par.epoch_ref)*par.day2tu;
	M0 = ephref0(6)+2*pi*epoch_diff0/par.T0;
	th0 = M2theta(M0, ephref0(2));
	kep0 = [ephref0(1:5), th0];
	rv0 = kep2cart(kep0,1.0);

	%destination
	epoch_diff1 = (mjd+tof-epoch_ref)*par.day2tu;
	M1 = ephref1(6)+2*pi*epoch_diff1/par.T1;
	
	th1 = M2theta(M1, ephref1(2));
	kep1 = [ephref1(1:5), th1];
	rv1 = kep2cart(kep1,1.0);
	toff = tof*par.day2tu;
	dvtot = [];
	for lw=0:1
        [vsc0,vsc1,a,p,theta,iter]=lambertI(rv0(1:3)',rv1(1:3)',toff,1.0,lw);
        dv0v = vsc0-rv0(4:6)';
        dv1v = vsc1-rv1(4:6)';
% 	        a1 = atan2d(norm(cross(dv0v,rv0(4:6)')),dot(dv0v,rv0(4:6)'))
% 	        a2 = atan2d(norm(cross(dv1v,rv0(4:6)')),dot(dv1v,rv1(4:6)'))
		dv0 = norm(dv0v)*par.au/par.tu;
		dv1 = norm(dv1v)*par.au/par.tu;
        dvtot(lw+1) = dv0+dv1;
        % dvtot = [dvtot,(dv0+dv1)];
 		% dv(end+1) = dvtot;
	end
	dv = min(dvtot); %km/s
	dvlog = log10(dv);
return