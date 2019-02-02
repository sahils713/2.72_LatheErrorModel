classdef lathe
    %LATHE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        E_st = 200E3; %[MPa] Modulus of Elasticity of Steel
        E_al = 70E3; %[MPa] Modulus of Elasticity of Aluminum
        alpha_st = 12E-6; %[/K] Linear thermal expansion coefficient of St
        alpha_al = 22E-6; %[/K] Linear thermal expansion coefficient of Al
    end
    
    properties % PART
        L1 = 40; %[mm] Length of the Part
        D1 = 12.5; %[mm] Diameter of the Part
        delT1 = 8; %[deg C]
        mat1 = 'Al'; % Material Selection
    end
    
    properties % CHUCK
        L2 = 30; %[mm] Length of Chuck
        delT2 = 8; %[deg C]
        al2 = 12E-6; %[/K] Linear thermal expansion coefficient of St;
    end

    properties % SPINDLE SHAFT
        D3 = 19; %[mm] Diameter of the spindle shaft
        L3 = 21; %[mm] Spindle shaft length from chuck to closest bearing
        delT3 = 10; %[deg C] Temperature Rise during Operation
        shaftEndTreatment = 'simple-support'; % Support Type ('simple-support' or 'cantilever')
        E3 = 200E3; %[MPa] Modulus of Elasticity of Steel
        al3 = 12E-6; %[/K] Linear thermal expansion coefficient of St
    end

    properties % SPINDLE SHAFT BEARINGS
        shaftBearingSpacing = 61; %[mm] Spacing between bearings
        radialBearingStiffness = 284400; %[N/mm] - Radial Bearing Stiffness
        axialBearingStiffness = 52540; %[N/mm] - Radial Bearing Stiffness
        spindleMisalignAngleX = deg2rad(-0.002); %[Rad]
        spindleMisalignAngleY = deg2rad(0.002); %[Rad]
        spindleBearingRunout = 0.0001*25.4; %[mm]
        spindleRunoutMisalignment = deg2rad(10); %[deg] Angular misalignment between high runout points
    end
    
    properties % BELT
        beltTension = 10; %[N] Belt Tension
        beltDistance = 20; %[mm] Belt Distance
    end

    properties % HEADSTOCK
        L4 = 60; % [mm] x-Distance from rails to spindle-bearing axis
        L5 = 100; % [mm] y-Distance from rails to bearing axis
        b5 = 15; % [mm] - Thickness of headstock
        w5 = 86; % [mm] - Width of headstock
        E5 = 70E3; %[MPa] Modulus of Elasticity of Aluminum
        headStockTorsionStiffness = 9.3612e+07; %[N mm/rad] based on FEA
    end
    
    properties % RAIL SYSTEM
        totalRailLength = 220; % [mm] Total Rail Length (z-direction)
        railSpacing = 115; %[mm] Rail Spacing (x-direction)
        railBearingSpacing = 60; % [mm] Bearing-to-Bearing Spacing (in z-direction)
        thirdBearingOffset = 10; % [mm] Bearing on 2nd rail offset from centerline (in z-direction)
        d6 = 3/4*25.4; % [mm] Rail diameter
        E6 = 200E3; %[MPa] Modulus of Elasticity of Steel
        railEndTreatment = 'pin-pin';
        railStraightness = (0.0002)/(12); % [mm/mm] Straightness {McMaster}
        railBowOrientation = pi()/2; %[rad] Plane of Maximum Bow oriented about z-axis
        rail2BowOrientation = pi()/2; %[rad] Plane of Maximum Bow oriented about z-axis (Rail with Flexure)
        railMisalignAngleX = (0.0005*25.4)/220; %[rad] % Axis misalignment about x-axis
        railMisalignAngleY = (-0.0005*25.4)/220; %[rad] % Axis misalignment about y-axis
        rail2MisalignAngleX = (0.0005*25.4)/220; %[rad] % Axis misalignment about x-axis (Rail with Flexure)
    end

    properties % Z-AXIS LEAD SCREW
        dlScrew = 20; %[mm] Diameter of the Leadscrew
        E_lScrew = 200E3; %[MPa] Modulus of Elasticity of Steel
    end

    properties % CARRIAGE
        L7 = 80; %[mm] Carriage Height (y-direction)
    end
    
    % CROSS-FEED ASSEMBLY
    properties
        % L8 is a dependent variable based on the cutting tool stick-out
        E8_lscrew = 200000; %[MPa] Modulus of Elasticity of Steel for Lead Screw
        flexureXstiffness= 4 * 480/1000; %[N/mm] based on Victor's Calculations
        xFeedStiffnessX = 3000; %[N/mm] based on measurements
        flexureYstiffness = 732; %[N/mm] based on Victor's Calculations
        flexureZstiffness = 24600; %[N/mm] based on Victor's Calculations
        XFeedLeadScrewDiameter = 0.25*25.4; %[mm] Diameter of Leadscrew
        flexureMountXSpacing = (3.165 + 2.705)*25.4; %[mm] Spacing of Flexure Mounts in X
        flexureMountZSpacing = 1.5*25.4;  % [mm] Spacing of Flexure Mounts in Z
        xFeedMisalignmentX = (0.1E-3)/3; % [rad] Misalignment aboux x-axis. Baseline - 0.1
        xFeedMisalignmentY = deg2rad(0.29); % [rad] Misalignment about y-axis. Baseline - 0.04 deg.
        xFeedMisalignmentZ = (2.5E-3)/0.580;% [rad] Misalignment aboux z-axis. Baseline - 0.04 deg.
        xFeedErrorM_thX = deg2rad(-0.0005); % [rad] Error Motion in thY through travel in X
        xFeedErrorM_thY = deg2rad(0.0005); % [rad] Error Motion in thZ through travel in X
        xFeedErrorM_Y = -0.010;  % [mm] Error Motion in Y through travel in X
        xFeedErrorM_Z = -0.010; % [mm] Error Motion in Z through travel in X        
    end
    
    % CUTTING TOOL POST
    properties
        % L9 is a dependent variable based carriage height and headstock
        % height
        L12 = 30; %[mm] Cutting tool offset from centerline of flexure in -z direction
    end

    properties % CUTTING TOOL
        L10 = 10; % [mm] Cutting Tool Stick-Out (x-direction)
        w10 = 6; % [mm] Cutting Tool Width (assume square)
        E10 = 200E3; %[MPa] Modulus of Elasticity of Steel
        al10 = 12E-6; %[/K] Linear thermal expansion coefficient of St
        delT10 = 10; % Temperature rise
        cuttingDepth = 0; % [mm] Depth of cut (radial)
    end
    
    properties(Dependent)
        E1; % Modulus of Elasticity of Part
        al1; % Linear Thermal Expansion Coefficient of Part
        distFirstBearing %[mm] Headstock to First Bearing on rail (z-direction)
        L6 % Modeled length of rail in HTM
        L8; %[mm] Effective flexure length in x-direction
        L9% [mm] Cutting Tool Post Height (y-direction)
        L11% [mm] Offset from Cutting Tool to Third Bearing(on 2nd Rail)
    end
    
    methods
        % Modulus of Elasticity of the Part
        function modulus = get.E1(obj)
            switch obj.mat1
                case 'Al'
                    modulus = obj.E_al;
                case 'St'
                    modulus = obj.E_st;
            end
        end
        
        % Linear Thermal Expansion Coefficient of Part
        function coeff = get.al1(obj)
            switch obj.mat1
                case 'Al'
                    coeff = obj.alpha_al;
                case 'St'
                    coeff = obj.alpha_st;
            end
        end
        
        % Distance from Headstock to First Bearing on rail (z-direction)
        function dis = get.distFirstBearing(obj)
            dis = obj.L1 + obj.L2 + obj.L12 - 0.5*obj.railBearingSpacing;
        end
        
        % Modeled length of rail in HTM
        function length = get.L6(obj) 
            length = obj.L1 + obj.L2 + obj.L3 + ...
                obj.L11 + obj.L12;
        end
        
        %[mm] Effective flexure length in x-direction
        function length = get.L8(obj)
            length = obj.L4 - obj.L10 - obj.D1/2 + obj.cuttingDepth;
        end
        
        % [mm] Cutting Tool Post Height in y-direction
        function length = get.L9(obj)
            length =  obj.L5-obj.L7;
        end
        
        % [mm] Offset from Cutting Tool to Third Bearing in z-direction
        function length = get.L11(obj)
            length = obj.thirdBearingOffset;
        end
    end
    
end

