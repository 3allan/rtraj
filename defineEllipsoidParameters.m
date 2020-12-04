function [GM, Req, Rpo, f, e, w] = defineEllipsoidParameters(EllipsoidModel)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Define the parameters composing the ellipsoid model of Earth.
% 
%  Inputs:
% 
%    EllipsoidModel - Name of model for the ellipsoidal Earth.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                       "GRS80" - Provides the Geodetic Reference System
%                                 (GRS) model of 1980. This model is used
%                                 with the International Terrestrial
%                                 Reference System (ITRS) and maintained by
%                                 the International Earth Rotation and
%                                 Reference System's (IERS) to establish an
%                                 accurate means of measuring positions
%                                 near Earth's surface.
% 
%                       "WGS84" - Provides the World Geodetic System (WGS)
%                                 model of 1984. This model is perhaps the
%                                 most widely used geodetic model of Earth.
%                                 It is compatible with the EGM2008 gravity
%                                 model (as well as the resulting geoid
%                                 parameters) and the WMM2020 magnetic
%                                 model without needing any geodetic
%                                 corrections/conversions between other
%                                 models. All GPS satellites and
%                                 coordinates coming from GPS applications
%                                 reference the WGS84 ellipsoid.
% 
%                     "IERS03" - Provides the International Earth Rotation
%                                and Reference System's (IERS) model of
%                                2003. This model is the current best
%                                estimate from the IERS of Earth's
%                                reference ellipsoid. Despite this
%                                ellipsoid being the best current estimate,
%                                common applications such as indicating
%                                geographic locations usually do not use
%                                this model since GPS use WGS84.
% 
%                     "CUSTOM" - Provides a set of user-defined ellipsoid
%                                parameters.
% 
%  Outputs:
% 
%                GM - Gravitational parameter of the ellipsoid model. The
%                     gravitational parameter is the standard product of
%                     Newton's gravitational constant (G) with the mass of
%                     Earth (M). Note that neither of G nor M are known
%                     very well independently due to how weak gravitational
%                     fields are and the inability to measure Earth's mass
%                     directly.Their product, however, constituting the 
%                     gravitational parameter GM is known to much greater 
%                     precision.
%                     Size: 1-by-1 (scalar)
%                     Units: m3/s2 (cubic meters per squared second)
% 
%               Req - Equatorial radius of the ellipsoid model. The
%                     equatorial radius is equivalent to the semimajor axis
%                     in the sense of a purely mathematical ellipsoid. In
%                     terms of the actual Earth, Req corresponds to an
%                     average radius from Earth's center at all longitudes
%                     around the equator, which is strongly constant in
%                     practice due to Earth's rigidity and slow precession
%                     and nutation of its spin vector.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               Rpo - Polar radius of the ellipsoid model. The polar radius
%                     is equivalent to the semiminor axis in the sense of a
%                     purely mathematical ellipsoid. In terms of the actual
%                     Earth, Rpo corresponds to the average distance from
%                     Earth's center to its two poles. These two distances
%                     are essentially the same for the same reasons of Req.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 f - The flattening of the ellipsoid model. The flattening
%                     describes how squashed (on its poles) the ellipsoid
%                     is from a perfect sphere. A flattening of 0 is
%                     equivalent to the statement that the ellipsoid is
%                     actually a perfect sphere. A flattening of 1 is
%                     equivalent to the statement that the ellipsoid is
%                     actually a flat (2D) circle (of radius Req). Thus,
%                     the flattening f satisfies 0 < f < 1 for every
%                     ellipsoid and is defined by f = (Req - Rpo) / Req.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%                 e - The eccentricity of the ellipsoid model. The
%                     eccentricity describes how far away the focal ellipse
%                     (or simply the 2 foci in an arbitrary cross-section of
%                     the (3D) ellipsoid) is away from its geometric
%                     center. An eccentricity of 0 is equivalent to the
%                     statement that the ellipsoid is actually a perfect
%                     sphere (such that the focal ellipse is a single point
%                     at the origin - a distance of 0 from the center). An
%                     eccentricity of 1 is equivalent to the statement that
%                     the ellipsoid is actually a flat (2D) circle (of
%                     radius Req). Thus, the eccentricity e also satisfies
%                     0 < e < 1 for every ellipsoid and is defined by
%                     e2 = 2f - f2.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%                 w - The rotation rate of the ellipsoid model. The
%                     rotation rate is a measure of how fast the model
%                     spins about its vertical axis piercing through the
%                     (north) pole. A positive angular rate indicates a
%                     westward-to-eastward rotational motion as experienced
%                     on the actual Earth. The actual Earth does not spin
%                     with an axis exactly through the North pole (it
%                     precesses and nutates in reality so that rotation is
%                     actually almost never directly through the North
%                     pole) nor does it spin with a constant angular rate,
%                     but such a model provides the means to approximate an
%                     Earth-centered inertial coordinate system without
%                     leading into overly-excessive details for very short
%                     flight times.
%                     Size: 1-by-1 (scalar)
%                     Units: rad/s (radians per second)
% 

% Check if the input is of type string (or character))
checkInput(EllipsoidModel)

% Convert to lowercase
EllipsoidModel = lower(EllipsoidModel);

% Provide the parameters according to EllipsoidModel
switch EllipsoidModel
    case "grs80"
        GM = 3.986004415e14; % Gravitational parameter [m3/s2]
        Req = 6378137.0; % Equatorial radius [m]
        f = 1/298.257222101; % Flattening []
        e = sqrt(2*f - f^2); % Eccentricity []
        Rpo = Req*(1 - f); % Polar radius [m]
        w = 7.292115855e-5; % Rotation rate [rad/s]
    case "wgs84"
        GM = 3.986004415e14; % Gravitational parameter [m3/s2]
        Req = 6378137.0; % Equatorial radius [m]
        f = 1/298.257223563; % Flattening []
        e = sqrt(2*f - f^2); % Eccentricity []
        Rpo = Req*(1 - f); % Polar radius [m]
        w = 7.292115855e-5; % Rotation rate [rad/s]
    case "iers03"
        GM = 3.986004415e14; % Gravitational parameter [m3/s2]
        Req = 6378136.6; % Equatorial radius [m]
        f = 1/298.256420000; % Flattening []
        e = sqrt(2*f - f^2); % Eccentricity []
        Rpo = Req*(1 - f); % Polar radius [m]
        w = 7.292115855e-5; % Rotation rate [rad/s]
    case "custom"
        GM = 3.986004415e14; % Gravitational parameter [m3/s2]
        Req = 6378136.3; % Equatorial radius [m]
        w = 7.292115855e-5; % Rotation rate [rad/s]
        e = 0.0818191908426215; % Ellipticity []
        f = 1 - sqrt(1 - e^2); % Flattening []
        Rpo = Req*(1 - f); % Polar radius [m]
    otherwise
        error("Please indicate a valid Earth ellipsoid model.")
end