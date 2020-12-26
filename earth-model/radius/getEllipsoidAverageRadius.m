function Ravg = getEllipsoidAverageRadius(a, b, methodOfAverage)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Calculate the average radius of the ellipsoid having semimajor axis a and
% semiminor axis b.
% 
%  Inputs:
% 
%                 a - Equatorial (semimajor) radius of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 b - Polar (semiminor) radius of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%   methodOfAverage - Name of method for calculating an average radius of
%                     the ellipsoidal Earth model of equatorial (semimajor)
%                     radius a and polar (semiminor) radius b.
%                     Size: 1-by-1 (string)
%                     Units: N/A (N/A)
% 
%                    Permissible options are:
%                          "2D" - Calculates the average radius of the
%                                 ellipsoid in the sense of the typical
%                                 integral-defined average. The integrand
%                                 is taken to be the radius r(t) of the 2D
%                                 vertical cross-section of the ellipsoid,
%                                 obtained by cutting the ellipsoid with
%                                 any vertical plane containing the z-axis.
%                                 The radius r(t) is therefore taken to be
%                                 r(t) = 1 / sqrt(1 - e2 cos(t)), where e2
%                                 is the ellipsoid's squared eccentricity 
%                                 and 0 < t < 2pi is the polar angle.
% 
%                          "++" - Calculates the average radius of the
%                                 ellipsoid in the sense of an algebraic
%                                 average of the equatorial
%                                 (semimajor) and polar (semiminor) radii
%                                 of the 2D ellipse formed in the vertical
%                                 cross-section of the 3D ellipsoid.
% 
%                        "RMS2" - Calculates the average radius of the
%                                 ellipsoid in the sense of a root-mean-
%                                 square (RMS) sum of the equatorial
%                                 (semimajor) and polar (semiminor) radii
%                                 of the 2D ellipse formed in the vertical
%                                 cross-section of the 3D ellipsoid.
% 
%                          "3D" - Calculates the average radius of the
%                                 ellipsoid in the sense of the typical
%                                 integral-defined average. The integrand
%                                 is taken to be the radius r(t) of the 3D
%                                 ellipsoid obtained in spherical
%                                 coordinates of variables (r, t, l), where
%                                 r is the radial coordinate (the radius),
%                                 t is the colatitude (zenith) angle, and 
%                                 l is the longitude (azimuth) angle. The
%                                 radius is obtained independent of l as
%                                 r(t) = a sqrt(1 - e2) ...
%                                        sqrt(2 / (2 - e2(1 - cos(2t)))),
%                                 where e2 is the ellipsoid's squared
%                                 eccentricity, a is the equatorial
%                                 (semimajor) radius, and 0 < t < pi.
% 
%                         "+++" - Calculates the average radius of the
%                                 ellipsoid in the sense of an algebraic
%                                 average of the equatorial (semimajor) and
%                                 polar (semiminor) radii in each of the
%                                 three principal directions.
% 
%                        "RMS3" - Calculates the average radius of the
%                                 ellipsoid in the sense of a root-mean-
%                                 square (RMS) sum of the equatorial
%                                 (semimajor) and polar (semiminor) radii
%                                 in each of the three principal directions.
% 
%  Outputs:
% 
%              Ravg - Average radius of the ellipsoid as determined by the
%                     invoked method.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 

% Check if the input is of type string (or character))
checkInput(methodOfAverage)

switch methodOfAverage
    case "2D"
        e2 = 1 - (b / a)^2;
        Ravg = (ellipke(e2) + ellipke(e2 / (e2 - 1)) / sqrt(1 - e2)) / pi;
    case "++"
        Ravg = (a + b) / 2;
    case "RMS2"
        Ravg = sqrt((a^2 + b^2) / 2);
    case "3D"
        e2 = 1 - (b / a)^2;
        Ravg = 2 * b * ellipke(e2) / pi;
    case "+++"
        Ravg = (2*a + b) / 3;
    case "RMS3"
        Ravg = sqrt((2*a^2 + b^2) / 3);
    otherwise
        error("Please indicate a valid method for obtaining the ellipsoid's average radius.")
end