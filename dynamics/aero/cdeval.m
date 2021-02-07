function CD = cdeval(M, a, V, pars, stage)
% 
% Arsam Jafaryzad (arsam@vt.edu) - Dec 23, 2020
% 
% Evaluate the aerodynamic charateristics of  
% the rocket with respect to the angle of attack,
% Reynolds number, and Mach number.
%   
%    Inputs:
% 
%                 V - Velocity of the rocket given during flight. It is in
%                     the ENV frame, a non-intertial frame.
%
%                 a - Sound speed at current altitude of flight.
%
%                 M - Mach number at current altitude of flight.
%                  
%              pars - Conglomeration of input arguments related to various
%                     rocket parameters. Rocket parameters required mainly 
%                     include those related to determining the drag
%                     coefficient.
%                     Size: my pp, so very big
%                     Units: au
% 
%    Outputs:
% 
%                CD - The drag coefficient of the vehicle
%                            
% 

% No checks

%% Symbols

% Symbols needed for calculations

%                 A - Area
% 
%               A_B - Base Area
% 
%              A_Bf - Base Area of One Fin
% 
%              A_BN - Nose Base Area  
% 
%              A_BS - Conical Shoulder or Boattail Base Area
% 
%                Ar - Reference Area
%
%                Af - Area of Fin
% 
%                AT - Planform Area of One Exposed Fin Panel
% 
%              A_WB - Total Body Wetted Area
% 
%             A_WBA - Wetted Area of Booster Nose Afterbody
% 
%              A_WN - Nose Wetted Area
% 
%             A_WNA - Wetted Area of Sustainer Nose and the Cylindrical
%                     Afterbody Between the Nose and the Shoulder. if there
%                     is no Shoulder, the Afterbody is taken as teh Remainder
%                     of the Sustainer
% 
%              A_WS - Conical Shoulder or Boattail Wetted Area
% 
%             A_WSA - Wetted Area of the Shoulder or Boattail and the
%                     Portion of the Sustainer Aft of the Shoulder or
%                     Boattail
% 
%                AR - Aspect Ratio of Exposed Fin
% 
%                 a - Ambient Speed of Sound
% 
%                B3 - Third Order Busseman Irreversibility Coefficient
% 
%                CD - Drag coefficient 
% 
%              CD_B - Base Drag Coefficient
% 
%             CD_BT - Tail Trailing Edge Drag Coefficient
% 
%              CD_f - Skin Friction Drag Coefficient
% 
%             CD_LT - Tail Leading Edge Drag Coefficient
% 
%              CD_P - Pressure Drag Coefficient
% 
%              CD_T - Total Drag Coeffcient
% 
%             CD_TT - Tail Thickness Drag Coefficient
% 
%              CD_w - Wave Drag Coefficient
% 
%                Cf - Incmpressible Skin Friction Coefficient
% 
%              Cf_c - Compressible Skin Friction Coefficient
% 
%                Cl - Rolling Moment Coefficient
% 
%              Cl_p - Roll Damping Moment Coefficient Derivative
% 
%              Cl_d - Roll Forcing Moment Coefficient Derivative
% 
%                Cm - Pitch Moment Coefficient
% 
%              Cm_q - Pitch Damping Moment Coefficient Derivative
% 
%                CN - Normal Force Coefficient
% 
%              CN_a - Normal Force Coefficient Derivative
% 
%                CP - Pressure Coefficient    
% 
%                 c - Chord Length
% 
%              c_MA - Mean Aerodynamic Chord
% 
%               c_r - Fin Root Chord
% 
%               c_t - Fin Tip Chord
% 
%                 D - Aerodynamic Drag Force
% 
%                 d - Diameter
% 
%               d_N - Diameter of Nose Base
% 
%               d_s - Diameter of Shoulder or Boattail Base
% 
%                 F - Force
% 
%               f_B - Total Body Fineness Ratio
% 
%                FD - Diederich's Planform Correlation Parameter for Fins
% 
%               f_e - Equivalent Fineness Ratio of a Truncated Cone
% 
%               f_N - Nose Fineness Ratio
% 
%              f_NA - Fineness Ratio of the Nose and Nose Afterbody
% 
%              f_SA - Fineness Ratio of the Shoulder or Boattail and
%                     Shoulder Afterbody
% 
%                 h - Fin Trailing Edge Thickness
% 
%              K, k - Inteference Factors and General Constants
% 
%                K1 - First Order Busseman Coefficient
% 
%                K2 - Second Order Buseeman Coefficient
% 
%                K3 - Third Order Busseman Coefficient
% 
%                 L - Length
% 
%                Lr - Reference Length
% 
%              l_AF - Location of Shoulder or Boattail Measured from teh
%                     Nose Tip
% 
%              l_LR - Fin Root Leading Edge Wedge Length
% 
%               l_N - Nose Length
% 
%               l_0 - Total Body Length
% 
%               l_s - Shoulder or Boattail Length
% 
%               l_T - Location of Fin Leading Edge Intersection with Body,
%                     Measured from Nose Tip
% 
%              l_TR - Fin Root Trailing Edge Wedge Length
% 
%              lbar - Roll Moment About Longitudinal Axis
% 
%                 M - Mach Number
% 
%                 m - Cotangent
% 
%              mbar - Pitch Moment About Center of Gravity
% 
%                 N - Number of Fins (3 or 4)
% 
%                 n - Normal Aerodynamic Force
% 
%                Pr - Prandtle Factor
% 
%                Pf - Prandtle Factor for Swept Fins
% 
%                 p - Roll Rate
% 
%                 q - Pitch Rate
% 
%              qbar - Dynamic Pressure
% 
%                Re - Reynolds Number
% 
%                Rs - Surface Roughness Height
% 
%                rL - Fin Leading Edge Radius
% 
%                rt - Body Radius at Tail
% 
%                 S - Exposed fin Smispan Measure From Root Chord
% 
%                 t - Maximum Fin Thickness
% 
%                 V - Velocity
% 
%               V_B - Body Volume
% 
%                 w - Downwash Velocity
% 
%              X_CG - Center of Gravity Location Measured from Nose Tip
% 
%              X_LT - Distance Between Fin Root Leading Edge and the Fin
%                     Tip Leading Edge Measured Parallel to the Root
% 
%              Xbar - Longitudinal Center of Pressure Location Measured
%                     from Nose Tip
% 
%              Ybar - Spanwise Center of Pressure Locaiton Measured for the
%                     Longitudinal Axis
% 
%              Y_MA - Spanwise Location of c_MA Measured from the Root Chord
% 
%              alfa - Angle of Attack
% 
%              beta - Supersonic sqrt(M^2-1), subsonic sqrt(1-M^2)
% 
%               G_c - Midchord Line Sweep Angle
% 
%               G_q - Quarter Chord Sweep Angle
% 
%               G_L - Leading Edge Sweep Angle
% 
%               G_T - Trailing Edge Sweep Angle
% 
%               G_1 - First Region Boundary Sweep Angle
% 
%               G_2 - Second Regiond Boundary Sweep Angle
% 
%                 g - Ratio of Specific Heats
% 
%             ksi_L - Fin Leading Edge Wedge Half Angle in a Plane Parallel
%                     to the Flow
% 
%             ksi_T - Fin Trailing Edge Wedge Half Angle in a Plane Parallel
%                     to the Flow
% 
%               eta - Incidence of Surface to Flow in Plane Parallel to the
%                     Flow
% 
%                 A - Fin Dihedral Angle
% 
%                 k - Fin Taper Ratio, c_t/c_r
% 
%                mu - Mach Angle
% 
%                nu - Ambient Kinematic Viscosity
% 
%               rho - Ambient Density
% 
%               sig - Nose Tip Half Angle
% 
%               tau - (s+r_t)/r_t ratio

% Subscripts of Symbols
% 
%                 B - Body
% 
%              B(T) - Body in Prsence of Tail
% 
%                cr - Critical Value 
% 
%                 N - Nose
% 
%                 r - Root
% 
%                 S - Conical Shoulder or Boattail
% 
%                 T - Tail
% 
%              T(B) - Tail In Presence of the Body
% 
%                 t - Tip
% 
%                 W - Wetted
% 
%                 l - Value for One Fin

%% Tail Normal Force Coefficient Derivative

beta = sqrt(1-M^2);

if M > 1
    beta = sqrt(M^2 - 1);
end

AR = 1;
G_c = 1;
Ar = 1;
k = 1;

CN_a0 = 2*pi/beta;
FD = (beta*AR)/cos(G_c);

% Normal force coefficient derivative of one fin
CN_a1 = (2*pi*AR*(Af/Ar))/(2+sqrt(4+((beta*AR)/(cos(G_c)))^2));

% Busemann's Third Order Expansion of Compressible Flow Equations
K1 = 2/beta;
K2 = ((k+1)*M^4-4*beta^2)/(4*beta^4);
K3 = (k+1)*M^8+(2*k^2-7*k-5)*M^6+10*(k+1)*M^4-12*M^2+8;
B3 = (k+1)*M^4*((5-3*k)*M^4+4*(k-3)*M^2+8)*(1/48);

A = 1;
CN_at = 2*CN_a1*cos(A)^2;

fincount = 4;
CN_aT = (fincount/2)*CN_a1;

%% Tail Center of Pressure Location


