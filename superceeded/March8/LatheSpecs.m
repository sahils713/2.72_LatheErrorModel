lathe = struct;

%% CONSTANTS
E_st = 200E3; %[MPa] Modulus of Elasticity of Steel
E_al = 70E3; %[MPa] Modulus of Elasticity of Aluminum
alpha_st = 12E-6; %[/K] Linear thermal expansion coefficient of St
alpha_al = 22E-6; %[/K] Linear thermal expansion coefficient of Al

%% PART INFORMATION
% PART
lathe.L1 = 50; %[mm] Length of the Part
lathe.D1 = 12.5; %[mm] Diameter of the Part
lathe.delT1 = 8; %[deg C]
lathe.E1 = E_al; % Modulus of Elasticity
lathe.al1 = alpha_al; % Linear Thermal Expansion Coefficient

% CHUCK
lathe.L2 = 30; %[mm] Length of Chuck
lathe.delT2 = 8; %[deg C]
lathe.al2 = alpha_st;

% SPINDLE SHAFT
lathe.D3 = 19; %[mm] Diameter of the spindle shaft
lathe.L3 = 20; %[mm] Spindle shaft length from chuck to closest bearing
lathe.delT3 = 5; %[deg C] Temperature Rise during Operation
lathe.al3 = alpha_st; % Linear Thermal Expansion Coefficient
lathe.shaftEndTreatment = 'simple-support'; % Support Type ('fixed-free' or 'cantilever')
lathe.E3 = E_st; % Stiffness of the Spindle Shaft


% SPINDLE SHAFT BEARINGS
lathe.shaftBearingSpacing = 70; %[mm] Spacing between bearings
lathe.bearingStiffness = 1000000; %[N/mm] - Bearing stiffness


% HEADSTOCK
lathe.L4 = 60; % [mm] x-Distance from rails to spindle-bearing axis
lathe.L5 = 100; % [mm] y-Distance from rails to bearing axis
lathe.b5 = 15; % [mm] - Thickness of headstock
lathe.w5 = 86; % [mm] - Width of headstock
lathe.E5 = E_al; % Modulus of Elasticity

% RAIL SYSTEM
lathe.totalRailLength = 220; % [mm] Total Rail Length (z-direction)
lathe.railSpacing = 115; %[mm] Rail Spacing (x-direction)
lathe.cutToolDistance = 50; % [mm] First Bearing (min z) to Cutting Tool (z-direction)
lathe.distFirstBearing = lathe.L1 + lathe.L2 -...
    lathe.cutToolDistance; %[mm] Headstock to First Bearing (z-direction)
lathe.railBearingSpacing = 83; % [mm] Bearing-to-Bearing Spacing (in z-direction)
lathe.thirdBearingOffset = 41.5; % [mm] 1st Bearing to only Bearing on 2nd rail (in z-direction)
lathe.d6 = 20; % [mm] Rail diameter
lathe.L6 = lathe.L1 + lathe.L2 + lathe.L3; % Modeled length in HTM
lathe.E6 = E_st; % Young's Modulus for Rail
lathe.railEndTreatment = 'fixed-fixed';
lathe.railStraightness = (3/16)/(10*12)/2; % [mm/mm] Straightness {McMaster}

% LEAD SCREW
lathe.dlScrew = 15; %[mm] Diameter of the Leadscrew
lathe.E_lScrew = E_st; % Young's Modulus for Leadscrew

% CARRIAGE
lathe.L7 = 80; %[mm] Carriage Height (y-direction)

% CUTTING TOOL POST
lathe.L9 = lathe.L5-lathe.L7; % [mm] Cutting Tool Post Height (y-direction)

% CUTTING TOOL
lathe.L10 = 20; % [mm] Cutting Tool Stick-Out (x-direction)
lathe.w10 = 6; % [mm] Cutting Tool Width (assume square)
lathe.E10 = E_st; % Young's Modulus of the Cutting Tool
lathe.al10 = alpha_st; % Linear Thermal Expansion Coefficient
lathe.delT10 = 0; % Temperature rise


% FLEXURE
lathe.L8 = lathe.L4-lathe.L10; % [mm] Flexure Length (x-direction)
lathe.flexureXstiffness= 4 * 480/1000; %[N/mm] based on Victor's Calculations
lathe.flexureYstiffness = 4 * 5.4E2; %[N/mm] based on Victor's Calculations
lathe.flexureZstiffness = 4 * 5.5E4; %[N/mm] based on Victor's Calculations