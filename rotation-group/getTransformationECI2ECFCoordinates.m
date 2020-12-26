function Tecf_eci = getTransformationECI2ECFCoordinates(t)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Obtain the linear transformation matrix T(t) converting components of
% vectors expressed in ECI coordinates to those expressed in ECF
% coordinates using the IAU-2000/2006 reduction method.
% 
%    Inputs:
% 
%                t  - Time
%                     Size: 1-by-1 (scalar)
%                     Units: JD UTC (Julian date in UTC)
% 
%    Outputs:
% 
%          Tecf_eci - Linear transformation matrix that transforms
%                     components of a vector expressed in the ECI frame to
%                     be expressed in the ECF frame by the linear
%                     transformation
%                                 x_ecf = Tecf_eci x_eci.
%                     Thus, the inverse transformation is obtained by
%                     a simple transpose operation (T is a member of SO3).
%                     Size: 3-by-3 (matrix - SO3)
%                     Units: - (unitless)
% 

% Convert this time t [JD] to UTC
tUTC = datetime(t, 'convertfrom', 'juliandate', 'timezone', 'utc');
tUTCvec = [tUTC.Year, tUTC.Month, tUTC.Day, tUTC.Hour, tUTC.Minute, tUTC.Second];
% Obtain transformation from ECI to ECF using the IAU-2000/2006 method and
% including 37 seconds as difference between TAI and UTC (10 + 27 leap
% seconds).
numberOfLeapSeconds = getLeapSeconds;
dTAIUTC = 10 + numberOfLeapSeconds;
Tecf_eci = dcmeci2ecef('IAU-2000/2006', tUTCvec, dTAIUTC);