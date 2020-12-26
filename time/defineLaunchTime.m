function [JDLaunch, tLaunchUTC, tLaunch] = defineLaunchTime(date, format, offset)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Define the initial time at which launch occurs by specifying a compatible
% pair (date, format) with a corresponding timezone offset indicating how 
% far ahead/behind the local area is from London, UK. Compatible here means 
% that the pair is accepted by datetime().
% 
%    Inputs:
% 
%              date - Specifies the initial time at which the vehicle
%                     launches. The format of it is a string/char spelling
%                     out the date according to `format'. Must be
%                     compatible with datetime().
%                     Size: 1-by-1 (string)
%                     Units: N/A (N/A)
% 
%            format - Accompanies `date' to tell the interpretter how to
%                     read the date correctly. Must be compatible with 
%                     datetime().
%                     Size: 1-by-1 (string)
%                     Units: N/A (N/A)
% 
%            offset - Accompanies `date' to tell the interpretter how to 
%                     read the time correctly. If no timezone is specified,
%                     then the timezone defaults to Coordinated Universal 
%                     Time (UTC). This time system has no hour offset from 
%                     UT1 or GMT, which are all observed in London, UK.
%                     Thus, specifying a time in `date' but no hour offset
%                     indicates that the time is interpreted to be local in
%                     London, UK. The offset must be in a form like
%                     '-04:00' to indicate a timezone 4 hours behind UTC or
%                     '+06:30' to indicate a timezone 6.5 hours ahead of
%                     UTC.
%                     Size: 1-by-1 (string)
%                     Units: N/A (N/A)
% 
%    Outputs:
% 
%          JDLaunch - Julian day time of launch according to the UTC time
%                     system. Note that Julian days count the days that
%                     have passed since Jan 1, 4713 BC 12:00:000 UT (noon
%                     and not midnight in London, UK).
%                     Size: 1-by-1 (scalar)
%                     Units: JD UTC (Julian days in UTC)
% 
%        tLaunchUTC - Time of launch according to the UTC time system.
%                     Size: 1-by-1 (datetime)
%                     Units: Time (UTC)
% 
%           tLaunch - Time of launch in the local timezone. Thus, tLaunch
%                     simply returns the indicated date, but as a datetime
%                     object.
%                     Size: 1-by-1 (datetime)
%                     Units: Time (local timezone)
% 

% Check that format and timezone are strings/chars
checkInput(date, format, offset);

% Obtain launch time in local time (same as specified date)
tLaunch = datetime(date, 'format', format, 'timezone', offset);

% Obtain launch time in UTC
tLaunchUTC = datetime(tLaunch, 'timezone', 'utc');

% Obtain Julian date of launch in UTC - does not account for any leap
% seconds.
JDLaunch = juliandate(tLaunchUTC);
