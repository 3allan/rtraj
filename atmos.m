function [d, p, T] = atmos(model, geomH, varargin)
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
%                                atmosphere model, is not precisely a
%                                standard atmosphere since the sphere is
%                                not the ellipsoid.
% 
%                        "ISA" - Provides the International Standard
%                                Atmosphere. This atmosphere is completely
%                                analogous to the 1976
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
%                                [1] http://braeunig.us/space/...
%                                    ...atmmodel.htm#refatmos
% 
%                 "MIL-STD210" - Provides a nonstandard atmosphere obtained
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
%                                [2] >> doc atmosnonstd
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
%                                successor to MIL-STD210.
%                                Reference:
%                                [3] >> doc atmosnonstd
%                                Aerospace toolbox only.
% 
%                   "CIRA1986" - Provides the corrected version of the 
%                                (C)OSPAR* (I)nternational (R)eference
%                                (A)atmosphere of 1986. This atmosphere
%                                model varies with both time (time-of-year)
%                                and location (latitude). Thus, conditions
%                                at a latitude are the same regardless of
%                                physical time (time of day). It is defined
%                                from 0 to 120,000 meters above mean sea
%                                level (MSL) and for geodetic latitudes 
%                                not exceeding +/-80 degrees.
%                                *Committee on Space Research (COSPAR)
%                                Reference:
%                                [4] >> doc atmoscira
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
%                                [5] >> doc atmosnrlmsise00
%                                [6] https://ccmc.gsfc.nasa.gov/...
%                                    ...modelweb/models/nrlmsise00.php
%                                [7] https://www.nrl.navy.mil/ssd/branches/...
%                                    ...7630/modeling-upper-atmosphere
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
%                                [8] https://sol.spacenvironment.net/jb2008/
%                                Requires gfortran.
% 
%                              * A (brief) list of other atmosphere (and
%                                other) models may be found here:
%                                [9] https://ccmc.gsfc.nasa.gov/modelweb/
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
%    Outputs:
% 
%                 d - Air/gas density.
%                     Size: n-by-1 (vector)
%                     Units: kg/m3 (kilogram per cubic meter)
% 
%                 p - Air/gas pressure.
%                     Size: n-by-1 (vector)
%                     Units: kg/m3 (kilogram per cubic meter)
% 
%                 T - Environmental temperature.
%                     Size: n-by-1 (vector)
%                     Units: kg/m3 (kilogram per cubic meter)
% 

% Check that model is of type string/char
checkInput(model)

% Uppercase
model = upper(model);

% Implement atmosphere models
switch model
    case "JET"
        % Define the Jet model algebraically by introducing standard US
        % 1976 values
        Rair = 287.058;
        [g0, Reff, p0, T0, S, Y] = getUS76BaseValues;
        
        % Don't bother converting to geopotential height - only 19 m
        % difference at 11 km altitude (~0.22% incurred error in doing so)
        geopH = geomH;
        
        % Calculate temperature according to linear variation with
        % geopotential altitude
        deltaT = dTdgeopH*geopH;
        T = T0 + deltaT;
        
        % Calculate pressure according to hydrostatic equation
        p = p0 * (1 + deltaT/T0)^()
        
        d0 = 1.225;  % Base
        
        