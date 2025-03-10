tic;

%% Script to mexify the spaceART toolbox
disp(' ');
disp(' ');
disp('-------------------------------------------------');
disp('Creation of the mex files of the spaceART_toolbox');
disp('-------------------------------------------------');

%% Extension and thus folder where the mex file shall be saved:
ext = mexext;
if (strfind(ext,'a')>0)
    dash = '/';
else
    dash = '\';
end

%% Reference frame conversions
fprintf('\n');
outputDir = 'conversion';
disp('REFERENCE FRAME CONVERSION FUNCTIONS');
fprintf('kep2cart: ');
keyboard
str = ['mex src' dash 'mexKep2cart.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('cart2kep: ');
str = ['mex src' dash 'mexCart2kep.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('car_bplT: ');
str = ['mex src' dash 'mexCar_bplT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('euler_axis_angle: ');
str = ['mex src' dash 'mexEuler_axis_angle.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('rth_carT: ');
str = ['mex src' dash 'mexRth_carT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('car_rthT: ');
str = ['mex src' dash 'mexCar_rthT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('tnh_carT: ');
str = ['mex src' dash 'mexTnh_carT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('car_tnhT: ');
str = ['mex src' dash 'mexCar_tnhT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('tnh_radecT: ');
str = ['mex src' dash 'mexTnh_radecT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('radec_tnhT: ');
str = ['mex src' dash 'mexRadec_tnhT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('car_radecT: ');
str = ['mex src' dash 'mexCar_radecT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('radec_carT: ');
str = ['mex src' dash 'mexRadec_carT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('rth_tnhT: ');
str = ['mex src' dash 'mexRth_tnhT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('tnh_rthT: ');
str = ['mex src' dash 'mexTnh_rthT.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Time conversions
fprintf('\n');
outputDir = ['conversion' dash 'time'];
disp('TIME CONVERSION FUNCTIONS');
fprintf('hms2fracday: ');
str = ['mex src' dash 'mexHms2fracday.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('fracday2hms: ');
str = ['mex src' dash 'mexFracday2hms.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('jd2mjd: ');
str = ['mex src' dash 'mexJd2mjd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd2jd: ');
str = ['mex src' dash 'mexMjd2jd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('jd2mjd2000: ');
str = ['mex src' dash 'mexJd2mjd2000.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd20002jd: ');
str = ['mex src' dash 'mexMjd20002jd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd2mjd2000: ');
str = ['mex src' dash 'mexMjd2mjd2000.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd20002mjd: ');
str = ['mex src' dash 'mexMjd20002mjd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('jd2date: ');
str = ['mex src' dash 'mexJd2date.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('date2jd: ');
str = ['mex src' dash 'mexDate2jd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd2date: ');
str = ['mex src' dash 'mexMjd2date.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('date2mjd: ');
str = ['mex src' dash 'mexDate2mjd.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('mjd20002date: ');
str = ['mex src' dash 'mexMjd20002date.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('date2mjd2000: ');
str = ['mex src' dash 'mexDate2mjd2000.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Ephemerides
fprintf('\n');
outputDir = 'ephemerides';
disp('EPHEMERIDES FUNCTIONS');
fprintf('EphSS: ');
str = ['mex src' dash 'mexEphSS.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('EphSS_kep: ');
str = ['mex src' dash 'mexEphSS_kep.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('moon_eph: ');
str = ['mex src' dash 'mexMoon_eph.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('NeoEphemeris: ');
str = ['mex src' dash 'mexNeoEphemeris.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('uplanet: ');
str = ['mex src' dash 'mexUplanet.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('astro_constants: ');
str = ['mex src' dash 'mexAstro_constants.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Dynamic equations
fprintf('\n');
outputDir = 'dynamics';
disp('DYNAMIC EQUATIONS FUNCTIONS');
fprintf('d2b_eq: ');
str = ['mex src' dash 'mexD2b_eq.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_m: ');
str = ['mex src' dash 'mexD2b_eq_m.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_mvarIsp: ');
str = ['mex src' dash 'mexD2b_eq_mvarIsp.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_Ttmod: ');
str = ['mex src' dash 'mexD2b_eq_Ttmod.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_u: ');
str = ['mex src' dash 'mexD2b_eq_u.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_u2: ');
str = ['mex src' dash 'mexD2b_eq_u2.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d2b_eq_utmod: ');
str = ['mex src' dash 'mexD2b_eq_utmod.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d3b_eq_PlanetSun_m: ');
str = ['mex src' dash 'mexD3b_eq_PlanetSun_m.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('d4b_eq_EarthSunMoon_wsb: ');
str = ['mex src' dash 'mexD4b_eq_EarthSunMoon_wsb.c src' dash 'dynamics.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Orbital transfers
fprintf('\n');
outputDir = 'orbital_transfers';
disp('ORBITAL TRANSFERS FUNCTIONS');
fprintf('dv_insertion: ');
str = ['mex src' dash 'mexDv_insertion.c src' dash 'orbital_transfers.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('semiperiod_hohm: ');
str = ['mex src' dash 'mexSemiperiod_hohm.c src' dash 'orbital_transfers.c src' dash 'ephemerides.c src' dash 'conversions.c src' dash 'astro_constants.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Relative motion
fprintf('\n');
outputDir = 'relative_motion';
disp('RELATIVE MOTION FUNCTION');
fprintf('lrom: ');
str = ['mex src' dash 'mexLrom.c src' dash 'relative_motion.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Swingby
fprintf('\n');
outputDir = 'swingby';
disp('SWINGBY FUNCTION');
fprintf('swingby: ');
str = ['mex src' dash 'mexSwingby.c src' dash 'swingby.c src' dash 'conversions.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%% Tools
fprintf('\n');
outputDir = 'tools';
disp('TOOLS FUNCTIONS');
fprintf('cartprod: ');
str = ['mex src' dash 'mexCartprod.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');
fprintf('qck: ');
str = ['mex src' dash 'mexQck.c src' dash 'math_utils.c -outdir ' outputDir ' -DMEXCOMPILE'];
eval(str)
fprintf('done\n');

%%
disp(' ');
disp(['All mex files created and saved in ' num2str(ceil(toc)) ' seconds.']);
disp('----------------------------------------------');
disp(' ');
disp(' ');
