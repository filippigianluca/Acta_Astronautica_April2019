% 
% 1. Open MSIS_atm_model.m and change the table of values, if required
% (data from http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi)
% 2. Open Chebyshev.m and change the name of the file where the chebyshev
% coefficients will be saved (save_filename)
% 3. Run Chebyshev.m
% 4. In Propagator_Earth.m change the name of file to be loaded for the
% definition of the coefficients (line ~373 in Propagator Earth.m)
% according to the name used in Chebyshev.m (save_filename)
% 5. Run Propagator_Earth.m