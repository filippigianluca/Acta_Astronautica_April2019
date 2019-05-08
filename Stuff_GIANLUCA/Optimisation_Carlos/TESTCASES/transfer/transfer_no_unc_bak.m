mu = 1.32712440018e11;
au = 1.49597870691e8;
tu = sqrt(au^3/mu);

tu2day = tu/86400;
day2tu = 1/tu2day;

epoch_ref = 58200;
% ephemeris [a e i W w M] in AU and degrees
ephref0 = [2.767046248500289	.07553461024389638	10.5935097971363	80.30991865594387   73.11534200131032	 352.2304611765882];
ephref1 = [2.772766124346009	.2305054995905846	34.83687567654144   173.0838008845752	310.0063356539786	 334.3231878791698];
% ephemeris adim
ephref0(3:6) = deg2rad(ephref0(3:6));
ephref1(3:6) = deg2rad(ephref1(3:6));

% periods
T0 = 2*pi*sqrt(ephref0(1)^3);
T1 = 2*pi*sqrt(ephref1(1)^3);


figure()
hold on
dv = [];
mjd0 = 58300+600;
tf0 = 0+200;
for mjdi = 1:1000
	mjd = mjd0+mjdi;
	epoch_diff0 = (mjd-epoch_ref)*day2tu;
	M0 = ephref0(6)+2*pi*epoch_diff0/T0;
	th0 = M2theta(M0, ephref0(2));
	kep0 = [ephref0(1:5), th0];
	rv0 = kep2cart(kep0,1.0);

	M1 = ephref0(6)+2*pi*epoch_diff0/T1;
	th1 = M2theta(M1, ephref1(2));
	kep1 = [ephref1(1:5), th1];
	rv1 = kep2cart(kep1,1.0);

	% plot3(rv0(1),rv0(2),rv0(3),'k.','MarkerSize',0.5);
	% plot3(rv1(1),rv1(2),rv1(3),'r.','MarkerSize',0.5);

	% if mod(mjdi,50)==0
% 		plot3(rv0(1),rv0(2),rv0(3),'k.','MarkerSize',10);
% 		plot3(rv1(1),rv1(2),rv1(3),'r.','MarkerSize',10);
% 		quiver3(rv0(1),rv0(2),rv0(3),rv0(4),rv0(5),rv0(6),'k');
% 		quiver3(rv1(1),rv1(2),rv1(3),rv1(4),rv1(5),rv1(6),'r');
% 		drawnow

	% end

	for tof = 1:200

		epoch_diff1 = (mjd+tof-epoch_ref)*day2tu;
		
		M1 = ephref1(6)+2*pi*epoch_diff1/T1;
		
		th1 = M2theta(M1, ephref1(2));
		kep1 = [ephref1(1:5), th1];
		rv1 = kep2cart(kep1,1.0);
		ERROR = 1;
		toff = tof*day2tu;

		lw=vett(rv0(1:3),rv1(1:3));
		lw=sign(lw(3));
		lw = (lw ~= 1);
		[vsc0,vsc1,a,p,theta,iter]=lambertI(rv0(1:3)',rv1(1:3)',toff,1.0,lw);

		% [A,P,E,ERROR,vsc0,vsc1,TPAR,THETA] = lambertMR(rv0(1:3),rv1(1:3),toff,1.0,0,0,0);
		% if ERROR~=0
		% 	ERROR
		% 	vsc0
		% 	vsc1
		% end

		% [vsc0, vsc1, extremal_distances, exitflag] = lambert(rv0(1:3), rv1(1:3), toff, 0, 1);

        dv0v = vsc0-rv0(4:6)';
        dv1v = vsc1-rv1(4:6)';
        a1 = atan2(norm(cross(dv0v,rv0(4:6)')),dot(dv0v,rv0(4:6)'))
        a2 = atan2(norm(cross(dv1v,rv0(4:6)')),dot(dv1v,rv1(4:6)'))
		dv0 = norm(dv0v)*au/1e3/tu;
		dv1 = norm(dv1v)*au/1e3/tu;
        dvtot = (dv0+dv1);
 		% dv(end+1) = dvtot;

        dv(mjdi,tof) = dvtot; %km/s
	end

% 	if mod(mjdi,50)==0
% 
% 		save('dv_transfer_5','dv');
% 
% 	end
end