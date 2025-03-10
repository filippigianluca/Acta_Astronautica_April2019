classdef CelestialBody 
    % CELESTIALBODY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        name % Name of the Celestial Body
        a    % Semimajor axis [AU]     
        e    % Eccentricity
        i    % Inclination [deg]     
        OM   % Asc. Node/raan [deg]
        W    % Arg. Perigee [deg]     
        M0   % Mean anomoly, M at time given t0 [deg]
        t0   % Time at which Mo is given [MJD2000]   
        mag  % Abs Magnitude (i.e. intrinsic brightness)             
        mass % Mass [kg]        
    end
    
    methods
        
        % Default CelestialBod Constructor
        %  
        % Params: 
        %  a   Semimajor axis [AU]     
        %  e   Eccentricity 
        %  i   Inclination [deg]     
        %  OM  Asc. Node/raan [deg]
        %  W   Arg. Perigee [deg]     
        %  M0  Mean anomoly, M at time given t0 [deg]
        %  t0  Time at which Mo is given [MJD2000]  
        function celestial_obj = CelestialBody(name,a,e,i,om,w,M0,t0) 
           celestial_obj.name = name;
           celestial_obj.a    = a; 
           celestial_obj.e    = e;
           celestial_obj.i    = i;
           celestial_obj.OM   = om;
           celestial_obj.W    = w;
           celestial_obj.M0   = M0;
           celestial_obj.t0   = t0;
        end
           
        
        % Get the keplerian elements 
        %
        % kep
        %   - a   Semimajor axis [AU]     
        %   - e   Eccentricity 
        %   - i   Inclination [deg]     
        %   - OM  Asc. Node/raan [deg]
        %   - W   Arg. Perigee [deg]     
        %   - M0  Mean anomoly, M at time given t0 [deg]
        %   - t0  Time at which Mo is given [MJD2000] 
        function kep = getKeplerianElements(obj)
            
            % Create a new Keplerian Element Object                 
            kep = KeplerianElements(obj.a,obj.e,obj.i,obj.OM, obj.W, obj.M0, obj.t0);
            
        end
                
    end
    
end

