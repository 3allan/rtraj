function leapSeconds = getLeapSeconds
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Get the net number of leap seconds inserted into UTC due to the
% nonuniformity of Earth's rotation causing UT1 and UTC to deviate from
% each other by more than about 0.6 seconds. Leap seconds can be either 
% positive (resulting in, for example, Dec 31, 20xx 12:59:60) or likewise
% negative. A negative leap second takes away from the total number of 
% leap seconds while a positive leap second of course adds to the total 
% number of leap seconds. All previous leap seconds have been positive, 
% as +1 second has been added to UTC every time it has been decided that a
% leap second will be added to the calendar.
% 
%    Inputs:
% 
%                   -
% 
%    Outputs:
% 
%       leapSeconds - Net number of leap seconds added to UTC since the
%                     realization of TAI. The first leap second was added
%                     on December 31, 1971 11:59:60 UTC. The International
%                     Earth Rotation and Reference Systems Service (IERS)
%                     decides on whether a leap second will be added into
%                     the current year every 6 months. It is only possible
%                     to add a leap second on the very last second of
%                     either June 30 or December 31 of any given year.

% Provide the number of leap seconds added to UTC.
leapSeconds = 27;