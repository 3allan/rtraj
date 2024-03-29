function [d, T, varargout] = atmos(model, geomH, varargin)
% 
% Matt Werner (m.werner@vt.edu) - Dec 5, 2020
% 
% Estimate the air/gas properties (density, pressure, and temperature) at a
% specified height above mean sea level or, more generally, position and
% time depending on the complexity of the model used. Earth's atmospheric/
% gas still exist even at very high altitudes (thousands of kilometers), so
% density, pressure, and (environmental) temperature are always defined
% eeven though the particle count per cubic meter in such orbital regions
% is low.
% 
%    Inputs:
% 
%             model - Indication of which atmosphere model to implement.
%                     Note that different models require different amounts
%                     of input.
%                     Size: 1-by-1 (string)
%                     Units: - (N/A)
% 
%                    Permissible options are:
%                        "Jet" - Provides the Jet Standard Atmosphere. This
%                                atmosphere model is defined from 0 to
%                                11,000 meters above mean sea level (MSL).
%                                The Jet atmosphere is independent of both
%                                time and longitude/latitude. Thus,
%                                conditions worldwide are the same
%                                regardless of physical time (time-of-day &
%                                time-of-year) or station location (equator
%                                or North pole). This model implicitly
%                                takes the Earth's shape to be a perfect
%                                sphere.
% 
%                       "US76" - Provides the World Meteorological 
%                                Organization (WMO)'s 1976 U.S. Standard
%                                Atmosphere. This atmosphere model is
%                                equivalent to the Jet Standard Atmosphere
%                                model for geometric altitudes below 11,000
%                                meters above mean sea level (MSL) and
%                                continues to provide standard estimates
%                                for density, pressure, and temperature for
%                                geometric altitudes below 86,000 meters
%                                above MSL. The 1976 U.S. Standard
%                                Atmosphere model is independent of both
%                                time and longitude/latitude. Thus,
%                                conditions worldwide are the same
%                                regardless of physical time (time-of-day
%                                & time-of-year) or station location
%                                (equator or North pole). This model
%                                implicitly takes the Earth's shape to be a
%                                perfect sphere.
% 
%                      "XUS76" - Provides the extended 1976 U.S. Standard
%                                Atmosphere. This atmosphere model is 
%                                equivalent to the WMO's 1976 U.S. Standard
%                                Atmosphere model for geometric altitudes 
%                                below 86,000 meters above mean sea level 
%                                (MSL) and continues to provide standard
%                                estimates for density, pressure, and 
%                                temperature for geometric altitudes below
%                                1,000 kilometers above mean sea level.
%                                The extended 1976 U.S. Standard Atmosphere 
%                                model is independent of both time and 
%                                longitude/latitude. Thus, conditions 
%                                worldwide are the same regardless of 
%                                physical time (time-of-day & time-of-year)
%                                or station (equator or North pole). This
%                                model implicitly takes the Earth's shape
%                                to be a perfect sphere.
% 
%                    "XUS76NG" - Provides the extended 1976 U.S. Standard
%                                Atmosphere, but modified as to let
%                                gravitational acceleration on the surface
%                                of the ellipsoid vary with geodetic
%                                latitude (normal gravity, where "normal"
%                                means that the direction of the gravity
%                                vector is perpendicular to the ellipsoid
%                                at all points on the ellipsoid). This
%                                atmosphere model is defined exactly the
%                                same (uses all the same coefficients) as
%                                the extended 1976 U.S. Standard Atmosphere
%                                model except that sea-level gravitational
%                                acceleration varies with geodetic
%                                latitude. This atmosphere model is
%                                independent of time, but does depend on
%                                position (geodetic latitude). Thus,
%                                conditions worldwide are the same
%                                regardless of physical time (time-of-day
%                                & time-of-year). This model uses varying
%                                gravity in latitude (ellipsoid), but uses
%                                the same data used to source the extended
%                                1976 U.S. Standard Atmosphere (XUS76)
%                                (sphere). The Earth models, and thus this
%                                atmosphere model, is not 'precisely' a
%                                standard atmosphere since the sphere is
%                                not the ellipsoid, but is close enough for
%                                all practical purposes.
%                                Reference:
%                                [01] NASA-TM-X-74335 (U.S. Standard
%                                     Atmosphere, 1976)
%                                     Doc ID: 19770009539
%                                [02] http://braeunig.us/space/atmmodel.htm#USSA1976
% 
%                        "ISA" - Provides the International Standard
%                                Atmosphere. This atmosphere is completely
%                                analogous to the 1976 U.S. Standard
%                                atmosphere (in that it is independent of
%                                time and longitude/latitude) except it
%                                serves as the international standard.
%                                [11] https://en.wikipedia.org/wiki/International_Standard_Atmosphere
% 
%         "Groves/Jacchia1971" - Provides a standard atmosphere that varies
%                                with both time (time-of-year) and position
%                                (altitude and latitude). This atmosphere
%                                model is defined from 0 to 2,000,000
%                                meters above mean sea level (MSL). The
%                                model is emperical by nature and assumes
%                                that atmospheric properties in the 
%                                southern hemisphere are offset by those in
%                                the northern hemisphere by 6 (six) months.
%                                It provides no information about the
%                                atmosphere for positions situated at or
%                                above/below +/-70 degrees geodetic
%                                latitude.
%                                Reference:
%                                [03] http://braeunig.us/space/atmmodel.htm#refatmos
% 
%                "MIL-STD210C" - Provides a nonstandard atmosphere obtained
%                                by the Department of Defense (DoD) in
%                                January 1987 from a portion of climatic 
%                                data of the worldwide air environment.
%                                This atmosphere model is defined from 0 to
%                                80,000 meters above mean sea level (MSL).
%                                The data used within this model do not
%                                include any locations above the Earth
%                                above/below +/-60 degrees geodetic
%                                latitude. The model is independent of time 
%                                and position (longitude/latitude). This 
%                                model implicitly takes the Earth's shape 
%                                to be a perfect sphere. This model is the
%                                predecessor to MIL-HDBK310.
%                                Reference:
%                                [04] >> doc atmosnonstd
%                                Aerospace toolbox only.
% 
%                "MIL-HDBK310" - Provides a nonstandard atmosphere obtained
%                                by the Department of Defense (DoD) in
%                                June 1997 from a portion of climatic 
%                                data of the worldwide air environment.
%                                This atmosphere model is defined from 0 to
%                                80,000 meters above mean sea level (MSL).
%                                The data used within this model do not
%                                include any locations above the Earth
%                                above/below +/-60 degrees geodetic
%                                latitude. The model is independent of time 
%                                and position (longitude/latitude). This 
%                                model implicitly takes the Earth's shape 
%                                to be a perfect sphere. This model is the
%                                successor to MIL-STD210C.
%                                Reference:
%                                [04] >> doc atmosnonstd
%                                Aerospace toolbox only.
% 
%                   "CIRA1986" - Provides the corrected version of the 
%                                (C)OSPAR* (I)nternational (R)eference
%                                (A)atmosphere of 1986. This atmosphere
%                                model varies with both time (time-of-year)
%                                and location (latitude). Thus, conditions
%                                at a latitude are the same regardless of
%                                physical time (time of day). It is defined
%                                from 20,000 to 120,000 meters above mean
%                                sea level (MSL) and for geodetic latitudes 
%                                not exceeding +/-80 degrees.
%                                *Committee on Space Research (COSPAR)
%                                Reference:
%                                [05] >> doc atmoscira
%                                Aerospace toolbox only.
% 
%                 "NRLMSISE00" - Provides the NRL*MSISE** atmosphere model
%                                of 2001. This model provides a full 
%                                description of Earth's atmosphere, varying
%                                with time (year, time-of-year, and time-
%                                of-day) and location (latitude and 
%                                longitude), as well as with F10.7 
%                                (solar radio flux) values and magnetic 
%                                indices (specified either by default or
%                                directly provided.) It is defined from 0
%                                to 1,000,000 meters and has additional 
%                                options/outputs to provide per-molecule 
%                                (He, O, N2, O2, Ar, H, N) densities for an
%                                accurate description of satellite drag,
%                                including anomalous oxygen if desired.
%                                *Naval Research Lab (NRL)
%                                **Mass Spectrometer and Incoherant Scatter
%                                   Radar Exosphere (MSISE)
%                                References:
%                                [06] >> doc atmosnrlmsise00
%                                [07] https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
%                                [08] https://www.nrl.navy.mil/ssd/branches/7630/modeling-upper-atmosphere
%                                Aerospace toolbox only.
%                                
%                     "JB2008" - Provides the Jacchia-Bowman Thermospheric
%                                Density Model of 2008. This model provides
%                                a full description of Earth's atmosphere,
%                                varying with time (year, time-of-year, and
%                                time-of-day) and location (latitude and
%                                longitude) using current F10.7 (solar
%                                radio flux) values and magnetic field
%                                indices.
%                                References:
%                                [09] https://sol.spacenvironment.net/jb2008/
%                                Requires gfortran.
% 
%                              * A (brief) list of more atmosphere (and
%                                other) models may be found here:
%                                [10] https://ccmc.gsfc.nasa.gov/modelweb/
% 
%              geomH - Geometric (physical) height above mean sea level
%                     (MSL) at which to evaluate air/gas properties. The
%                     geopotential height (hgp) is related to the
%                     geometric, or physical, height h by
%                                 hgp = (Re h) / (Re + h),
%                     which does assume a spherical (as opposed to an
%                     ellipsoidal) Earth. This expression is used even in
%                     the more advanced atmosphere models, where Re is the
%                     MEAN radius of the Earth (Re = 6,371,010 meters).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%          varargin - See individual models for details.
%                     Size: ?
%                     Units: ?
% 
%    Outputs:
% 
%                 d - Air/gas density.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilogram per cubic meter)
% 
%                 T - Environmental temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%         varargout - Various outputs according to the model used. Each
%                     model assigns quantities into varargout at the end of
%                     their switch-case structure. In every case, the first
%                     argument of varargout (the third output) will ALWAYS
%                     be pressure with units of Pa (Pascals).
%                     Size: ?
%                     Units: ?
% 

% Check that model is of type string/char
checkInput(model)

% Uppercase
modelUpper = upper(model);

% Implement atmosphere models
switch modelUpper
    case "JET"
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for Jet model.")
        end
        % Enforce that the requested altitude is appropriate for the Jet
        % model
        if (geomH > 11019)
            error("Jet model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Don't bother converting to geopotential height - only 19 m
        % difference at 11 km altitude (~0.22% incurred error in doing so)
        Reff = 6356766; % Effective radius of Earth ([01] pg. 4) [m]
        geopH = geomH;
        
        % Define the Jet model algebraically by introducing standard US
        % 1976 values
        [d, T, p, c, mu] = atmosXUS76(geopH, 9.80665, Reff);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "US76"
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for US76 model.")
        end
        % Enforce that the requested altitude is appropriate for the US76
        % model
        if (geomH > 86000)
            error("1976 U.S. Standard Atmosphere model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth ([01] pg. 4) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call upon the 1976 U.S. Standard Atmosphere
        [d, T, p, c, mu] = atmosXUS76(geopH, 9.80665, Reff);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "XUS76"
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for XUS76 model.")
        end
        % Enforce that the requested altitude is appropriate for the XUS76
        % model
        if (geomH > 1e6)
            error("Extended 1976 U.S. Standard Atmosphere model is invalid for requested altitude (%1.0f m).", geomH)
        end
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth ([01] pg. 4) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call upon the 1976 U.S. Standard Atmosphere
        [d, T, p, c, mu] = atmosXUS76(geopH, 9.80665);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "XUS76NG"
        % Distribute varargin
        geodeticLatitude = varargin{1};
        Req = varargin{2};
        f = varargin{3};
        
        % Extended 1976 U.S. Standard Atmosphere with variations of
        % gravitational acceleration with respect to latitude (normal
        % gravity)
        [d, T, p, c, mu] = atmosXUS76NG(geomH, geodeticLatitude, Req, f);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "ISA"
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for ISA model.")
        end
        % Enforce that the requested altitude is appropriate for the US76
        % model
        if (geomH > 86000)
            error("International Standard Atmosphere model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth (same radius as US76) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call upon the International Standard Atmosphere (ISA)
        [d, T, p, c, mu] = atmosISA(geopH, gSL);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "Groves/Jacchia1971"
        error("Groves-Jacchia 1971 model not yet implemented.")
    case "MIL-STD210C"
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for MIL-STD210C model.")
        end
        % Enforce that the requested altitude is appropriate for the
        % MIL-STD210C model
        if (geomH > 2000000)
            error("MIL-STD210C model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth (same radius as US76) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call on the MIL-STD210C atmosphere model
        [T, c, p, d] = atmosnonstd(geopH, 'Profile', 'High density', ...
                                    '5%', 5, 'warning', '210c');
        
        % Calculate dynamic viscosity according to Sutherland's law
        mu = computeDynamicViscosity(T);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "MIL-HDBK310"
        % Identical segment to MIL-STD210C except that the specification in
        % atmosnonstd is '310' instead of '210c'
        
        % Enforce that at most 5 outputs are requested
        if (nargout > 5)
            error("Must request no more than 5 outputs (d, T, p, c, mu) for MIL-HDBK310 model.")
        end
        % Enforce that the requested altitude is appropriate for the
        % MIL-HDBK310 model
        if (geomH > 2000000)
            error("MIL-HDBK310 model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth (same radius as US76) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call on the MIL-STD210C atmosphere model
        [T, c, p, d] = atmosnonstd(geopH, 'Profile', 'High density', '10%', 5, 'warning', '310');
        
        % Calculate dynamic viscosity according to Sutherland's law
        mu = computeDynamicViscosity(T);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
    case "CIRA1986"
        % Enforce that at most 6 outputs are requested
        if (nargout > 6)
            error("Must request no more than 6 outputs (d, T, p, c, mu, wind) for CIRA1986 model.")
        end
        % Enforce that the requested altitude is appropriate for the
        % CIRA1986 model
        if (geomH > 120000)
            error("CIRA1986 model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Distribute varargin
        geodeticLat = varargin{1}; % Geodetic latitude [rad]
        month = varargin{2}; % Numeric representation of the month 
        %                           (Jan = 1, ..., Dec = 12)
        
        % Convert geodetic latitude to degrees
        geodeticLatdeg = rad2deg(geodeticLat);
        
        % Calculate geopotential altitude
        Reff = 6356766; % Effective radius of Earth (same radius as US76) [m]
        geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
        
        % Call on the CIRA 1986 model - Note that its direct outputs are
        % temperature, pressure and zonal wind (i.e. not density)
        [T, p, wind] = atmoscira(geodeticLatdeg, 'GPHeight', geopH, ...
                                  'monthly', month, 'warning');

        % Compute density, speed of sound, and dynamic viscosity on our own
        % using a constant specific gas constant of 287 J/kgK and the ideal
        % gas law
        Rair = 287;
        d = p / (Rair * T); % Density [kg/m3]
        c = sqrt(1.4*Rair*T); % Speed of sound [m/s]
        mu = computeDynamicViscosity(T);
        
        % Assign variable outputs
        varargout{1} = p; % Pressure [Pa]
        varargout{2} = c; % Speed of sound [m/s]
        varargout{3} = mu; % Dynamic viscosity [Pa*s]
        varargout{4} = wind; % Zonal wind [m/s]
    case "NRLMSISE00"
        % Enforce that at most 2 outputs are requested
        if (nargout > 9)
            error("Must request no more than 9 outputs (see atmos.m) for NRLMSISE00 model.")
        end
        % Enforce that the requested altitude is appropriate for the
        % NRLMSISE00 model
        if (geomH > 1000000)
            error("NRLMSISE00 model is invalid for requested altitude (%1.0f m).", geomH)
        end
        
        % Distribute varargin
        geodeticLat = varargin{1}; % Geodetic latitude [rad]
        longitude = varargin{2}; % Longitude [rad]
        year = varargin{3}; % Year
        dayOfYear = varargin{4}; % Day of year (Jan 01 = 1, Dec 31 = 365)
        UTseconds = varargin{5}; % Seconds passed since 00:00:000 UT in UK
        
        % Convert geodetic latitude and longitude to degrees
        geodeticLatdeg = rad2deg(geodeticLat);
        longitudedeg = rad2deg(longitude);
        
        % Define approximate local apparent solar time (NRLMSISE00
        % indicates that a more accurate expression may be used, but it is
        % of little importance/significance.
        localApparentSolarTime = UTseconds/3600 + longitudedeg/15;
        
        % Call on the NRLMSISE-00 atmosphere model 
        %   - No flags set (23 options available)
        %   - Default units output in the density slot is kg/m3
        [TxAlt, densities] = atmosnrlmsise00(geomH, geodeticLatdeg, ...
                               longitudedeg, year, dayOfYear, UTseconds, ...
                               localApparentSolarTime, 'Oxygen');
                                
        % Distribute the results
        T = TxAlt(:, 2); % Environmental temperature (temperature at altitude) [K]
        d = densities(:, 6); % Total mass density (for drag) [kg/m3]
        
        % Distribute the exospheric temperature and mass densities also
        varargout{1} = TxAlt(:, 1); % Exospheric temperature [K]
        varargout{2} = densities(:, 1); % Helium (He) density [1/m3]
        varargout{3} = densities(:, 2); % Oxygen (O) density [1/m3]
        varargout{4} = densities(:, 3); % Nitrogen (N2) density [1/m3]
        varargout{5} = densities(:, 4); % Oxygen (O2) density [1/m3]
        varargout{6} = densities(:, 5); % Argon (Ar) density [1/m3]
        varargout{7} = densities(:, 7); % Helium (H) density [1/m3]
        varargout{8} = densities(:, 8); % Nitrogen (N) density [1/m3]
        varargout{9} = densities(:, 9); % Anom. Oxygen density [1/m3]
    case "JB2008"
            error("Atmosphere not yet implemented - missing gfortran")
    otherwise
        error("Atmosphere model '%s' not found.", model)
end
        