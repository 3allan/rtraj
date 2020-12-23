# olvt

===============================================================================
                       The (R)ocket (Traj)ectory Program

                                                        ,:
                                                      ,' |
                                                     /   :
                                                  --'   /
                                                  \/ />/
                                                  / <//_\
                                               __/   /
                                               )'-. /
                                               ./  :\
                                                /.' '
                                              '/'
                                              +
                                             '
                                           `.
                                       .-"-
                                      (    |
                                   . .-'  '.
                                  ( (.   )8:
                              .'    / (_  )
                               _. :(.   )8P  `
                           .  (  `-' (  `.   .
                            .  :  (   .a8a)
                           /_`( "a `a. )"'
                       (  (/  .  ' )=='
                      (   (    )  .8"   +
                        (`'8a.( _(   (
                     ..-. `8P    ) `  )  +
                   -'   (      -ab:  )
                 '    _  `    (8P"Ya
               _(    (    )b  -`.  ) +
              ( 8)  ( _.aP" _a   \( \   *
            +  )/    (8P   (88    )  )
               (a:f   "     `"       `
			 	    rtraj
===============================================================================

rtraj is a Matlab program written to provide another means to simulate the 
trajectory of a rocket launching from Earth's surface. At its most fundamental
level, rtraj provides a simulation of the position and velocity of the 
rocket's mass center as measured relative to local level-vertical (LLV), also
called East-North-Vertical (ENV) or East-North-Up (ENU), coordinates. This
coordinate system is classified to be of the noninertial type and is treated
as such, which does not affect low-altitude (model) rockets significantly,
but becomes increasingly important with altitude and flight time. Beyond this
basic capability, rtraj provides a full 6 degree of freedom (6dof) description
of the flight vehicle to include body-orientation in addition to body position.

Running this version of rtraj requires ownership of several Matlab toolboxes 
and add-ons above the base product of Matlab. These additional packages are 
listed below. 
% ------------------------------------------------------------
@Toolboxes:
1. Aerospace toolbox
2. ...
% ------------------------------------------------------------
@Add-ons:
1. Ephemeris Data for Aerospace Toolbox (version 20.2.0)
% ------------------------------------------------------------
Without all of the above toolboxes and add-ons properly installed on the 
system, an error will simply ensue. rtraj will be unable to produce any 
results until all toolboxes and add-ons are added to the system.

rtraj is an open-source project and, as such, is required to be edited in the
source code to define the necessary parameters for a particular rocket until a 
GUI is implemented. It accepts multi-staged rockets; a single stage rocket has
parameters that are nearly all scalar-valued, but doubly-staged (and so on)
rockets have inputs that are nearly all vector/cell-valued to reflect upon the 
fact that the rocket contains multiple stages. All provided parameters and
outputs are stated in SI (MKS/metric) units.

The current implementation of rtraj does not yet include a
cache system and has not been optimized, so loading the various models may 
take a significant amount of time and computations may be fast or slow 
depending on the computer's amount of available RAM and computational 
power/ability. rtraj also fully relies upon being fed various models for
aerodynamic drag (drag coefficient) and propulsive thrust (thrust profile). It
may also be fed information regarding the mass flow rate, which is assumed to
take a linear variation if no mass flow rate profile is provided, and
information regarding the chamber pressure inside of the nozzle, which is not
required. Providing the chamber pressure profile yields extra information
about the gas flow in the nozzle.

% ------------------------------------------------------------
ASSUMPTIONS (not yet an exhaustive list)
1. ...
% ------------------------------------------------------------

Matt Werner (m.werner@vt.edu) - Dec 22, 2020