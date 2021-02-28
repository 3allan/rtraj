function [value, isterminal, direction] = odevents(t, x, earth, rocket, flags, s)
% 
% Matt Werner (m.werner@vt.edu) - Feb 12, 2021
% 
% Mark special occasions throughout the flight profile as an event. This
% function is called by the Matlab ODE solver upon specifying the event
% function in the ODE settings.
% 
%    Inputs:
% 
%                 t - Integration time.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%                 x - Integration state.
%                     Size: n-by-1 (vector)
%                     Units: ? (SI)
% 
%                 s - The current rocket's stage.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%             earth - Structure containing all information related to the
%                     Earth.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 
%            rocket - Structure containing all information related to the
%                     flight vehicle's intended flight profile.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 
%    Outputs:
% 
%             value - Determines whether an event is due according to its
%                     definition, which is (in general) dependent upon the
%                     earth, the rocket (its intended flight profile), the
%                     state, and time. 
%                     Note: An event occurs if any element crossing zero
%                           does so in the direction indicated by the
%                           corresponding element in 'direction'.
%                     Size: m-by-1 (vector)
%                     Units: ? (SI)
% 
%        isterminal - Decides whether or not the event should stop
%                     integration. Integration should be stopped at the end
%                     of an interval of ODE uniqueness. Such an end may
%                     correspond with any sudden discontinuity in the
%                     underlying dynamics determining the motion of the
%                     flight vehicle. For example, integration should be
%                     stopped upon staging the launch vehicle and at
%                     touchdown.
%                     Note: Stopping integration in the middle of an
%                           interval of uniqueness and continuing
%                           integration in the same interval does not
%                           affect uniqueness of that solution. Therefore,
%                           integration may be safely stopped and continued
%                           if uncertain that an action should halt
%                           integration.
%                     Size: m-by-1 (vector)
%                     Units: - (unitless)
% 
%         direction - Decides if a zero of 'value' is valid for marking an
%                     event. The zero's validity is determined by ensuring
%                     that it crosses zero in a particular specified
%                     direction. If the value crosses a zero but in the
%                     wrong direction, then there is no marking for an
%                     event. A direction of '0' marks events regardless of
%                     direction.
%                     Note: Scalars can cross 0 in two directions.
%                           1. (+1) - to + (upwards)
%                           2. (-1) + to - (downwards)
% 

% Make an event for fully leaving the launch rail and don't stop
% integration.
% Define the value as the difference between the position of the
% vehicle's base as measured in the rail frame and the launch rail
% length.
% As such, the direction must be upwards.
%
% Default value since this event only happens once
value(1) = -1;
isterminal(1) = 0;
direction(1) = 1;
% Ensure that this event may only be marked once - give it a maximum of 15
% seconds to clear the rail (Apollo 11 took 12 seconds)
if (t < 15)
    pos_rail = earth.T.rail_env*x(1:3);
    value(1) = pos_rail(3,1) - earth.launchsite.towerLength(1,1);
    isterminal(1) = 0;
    direction(1) = 1;
end

% Make an event for burnout and don't stop integration.
% Define the value as the difference between the current flight time and
% the current stage's burntime on top of the firing retarded time.
% As such, the direction must be upwards.
value(2) = t - ...
    rocket.motor.retardedTimes(s,1) - ...
    rocket.motor.burnTime(s,1);
isterminal(2) = 0;
direction(2) = 1;

% Make an event for staging and make it stop integration.
% Define the value as the difference between the current flight time and
% the last burnout time on top of the separation time.
% As such, the direction must be upwards.
value(3) = value(2) - rocket.delays.separation(s,1);
isterminal(3) = 1;
direction(3) = 1;