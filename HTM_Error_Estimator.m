function [delP, errorContributions, deflectionMatrix,HTMs, HTMs_d] = HTM_Error_Estimator(Fn, Ft, Fa, lathe,plotControl )
%% WORLD'S BEST LATHE HTM Error Estimator

% INPUTS
% Cutting Forces Fn, Fa, Fz in x, y, and -z direction [N]
% Lathe structure with dimensions, properties, and temperature details

% OUTPUTS
% Error P

% VERSION CONTROL
% Rev 3 (created - April 1 2016 - Sahil S)
% - Added carriage error motions into analysis

% Rev 2 (Created - Mar 21 2016 - Sahil S)
% - Added flexure and spindle misalignment

% Rev 1 (Created - Mar 19 2016 - Sahil S)
% - Added cutting-tool offset from center of flexure

% Rev 0 (Created - Mar 18 2017 - Sahil S)
% - Includes bearing stiffness from manufacturer
% - Added effect of belt tension on deflection
% - Rail bowing and misalignment included in error estimate

% Rev A (Created - Mar 04 2017 - Sahil S)

%% INITIALIZE HTM OUTPUTS
HTMs = struct; % Undeformed HTMs
HTMs_d = struct; % Deformed HTMs
HTMs_m = struct; % Magnified Error (for Visualization)


%% (A)PART
% Assumptions:
% - Moment due to off-axis axial load has not been included
% - Torsion about z-axis is a DOF
% - Part is cantilevered

% DIMS, PROPERTIES, and TEMP
[A1,I1] = shapeProps('circle',lathe.D1,[]); % Obtain Area and Inertia

% DEFLECTION
[dx_A,thy_A] = cantilever_F(Fn, lathe.L1, lathe.E1, I1);
[dy_A,thx] = cantilever_F(Ft, lathe.L1, lathe.E1, I1);
thx_A = -thx;
dz_A = -1*axial(Fa, lathe.L1, lathe.E1, A1) +...
    thermalExpansion(lathe.al1, lathe.L1, lathe.delT1); 
thz_A = 0;

% CONSTRUCT HTM
[HTMs.A, HTMs_d.A, HTMs_m.A] = HTM(-lathe.D1/2+lathe.cuttingDepth, 0, lathe.L1, dx_A, dy_A, dz_A, thx_A, thy_A, thz_A);

%% (B) CHUCK
% Assumptions:
% - Very stiff compared to other elements (NEED TO VERIFY)

% Calculate deflections
dz_B = thermalExpansion(lathe.al2, lathe.L2, lathe.delT2);
dx_B = 0; dy_B = 0; thx_B = 0; thy_B = 0; thz_B = 0;

% Construct non-deformed and deformed HTMs.
[HTMs.B, HTMs_d.B, HTMs_m.B] = HTM(0, 0, lathe.L2, dx_B, dy_B, dz_B, thx_B, thy_B, thz_B);

%% (C) SPINDLE SHAFT
% Assumptions:
% - Bearings can be treated as fixed-free support OR cantilever. Both cases
% investigated here.
% - Shaft deflections based on minimum diameter (conservative)

% DIMS, PROPERTIES, and TEMP
b = lathe.shaftBearingSpacing; %[mm] Spacing between bearings
[A3,I3] = shapeProps('circle',lathe.D3,[]); % Obtain Area and Inertia

% MOMENTS ON THE PART
My3  = Fn*(lathe.L1 + lathe.L2);
Mx3 = Ft*(lathe.L1 + lathe.L2);

% DEFLECTION
% Axial Loads + Thermal
dz_C = -1*axial(Fa, lathe.L3, lathe.E3, A3) +...
    thermalExpansion(lathe.al3, lathe.L3, lathe.delT3);
thz_C = 0;

switch lathe.shaftEndTreatment
    case 'cantilever'
        % From Point-Loads
        [dx_F,thy_F] = cantilever_F(Fn, lathe.L3, lathe.E3, I3);
        [dy_F,thx_F] = cantilever_F(Ft, lathe.L3, lathe.E3, I3);
        thx_F = -thx_F;

        % From Moment
        [dx_M,thy_M] = cantilever_M(My3, lathe.L3, lathe.E3, I3);
        [dy_M,thx_M] = cantilever_M(Mx3, lathe.L3, lathe.E3, I3);
        thx_M = - thx_M;
        
        % Sum directly for cantilever
        dx_C = dx_F + dx_M;       dy_C = dy_F + dy_M;
        thx_C = thx_F + thx_M;    thy_C = thy_F + thy_M;
        
    case 'simple-support'
        % Multiply by scale factors (see supporting calculations in binder)
        [Rb1x,Rb2x,dx_C,thy_C] = spindleDeflection(lathe.beltDistance,lathe.shaftBearingSpacing,lathe.L3, Fn, My3, 0, lathe.E3, I3);
        [Rb1y,Rb2y,dy_C,thx] = spindleDeflection(lathe.beltDistance,lathe.shaftBearingSpacing,lathe.L3, Ft, Mx3, lathe.beltTension, lathe.E3, I3);
        thx_C  = -thx;
end

% CONSTRUCT HTMS
[HTMs.C, HTMs_d.C, HTMs_m.C] = HTM(0, 0, lathe.L3, dx_C, dy_C, dz_C, thx_C, thy_C, thz_C);

%% (D) BEARING STIFFNESS
% Assumptions
% - Bearing stiffness to be equal for both bearings
% - Bearing stiffness assumed to be constant (VERIFY)

% DIMS AND PROPERTIES
K_b = lathe.radialBearingStiffness; %[N/mm] - Bearing stiffness
K_a = lathe.axialBearingStiffness; %[N/mm] - Bearing stiffness

% REACTION LOADS ON BEARINGS
% Calculated in previous section

% ANGULAR ERROR DUE TO BEARING RUNOUT (See derivation in Sahil's Notebook -
% Mar 21)
runoutAngularError = lathe.spindleBearingRunout/(2*b)*...
    (sin(deg2rad(45) + lathe.spindleRunoutMisalignment) - sin(deg2rad(45)));

% DEFLECTIONS AND ERRORS
dx_D = (Rb2x-Rb1x)/(2*K_b) + lathe.spindleBearingRunout; % Deflection in x-direction of bearing next to chuck
dy_D = (Rb2y-Rb1y)/(2*K_b) + lathe.spindleBearingRunout; % Deflection in y-direction of bearing next to chuck
dz_D = -Fa/K_a;
thx_D = -(Rb2y+Rb1y)/(K_b*b) + lathe.spindleMisalignAngleX - runoutAngularError; % Angular deflection about x-axis
thy_D = (Rb2x+Rb1x)/(K_b*b) + lathe.spindleMisalignAngleY + runoutAngularError; % Angular deflection about y-axis
thz_D = 0;

% CONSTRUCT HTMS
[HTMs.D, HTMs_d.D, HTMs_m.D] = HTM(0, 0, 0, dx_D, dy_D, dz_D, thx_D, thy_D, thz_D);

%% (E) HEADSTOCK
% Assumptions
% - Bolt Stiffnessess from bearing mount to headstock not included

% DIMS AND PROPERTIES
[~,Izz5] = shapeProps('rectangle',lathe.w5,lathe.b5); % Obtain Area and Inertia
[A5,Ixx5] = shapeProps('rectangle',lathe.b5,lathe.w5); % Obtain Area and Inertia

% MOMENTS
My5 = Fn*(lathe.L1 + lathe.L2 + lathe.L3);
Mx5 = Ft*(lathe.L1 + lathe.L2 + lathe.L3);

% DEFLECTIONS
[dx_E,thz] = cantilever_F(Fn, lathe.L5, lathe.E5, Izz5);
thz_E = -thz;
dy_E = axial(Ft, lathe.L4, A5, lathe.E5);

% Deflections due to Forces
[dz_F,thx_F] = cantilever_F(Fa, lathe.L5, lathe.E5, Ixx5);
thx_F = -thx_F; dz_F = -dz_F;
% Deflections due to Moments
[dz_M,thx_M] = cantilever_M(Mx5, lathe.L5, lathe.E5, Ixx5);
dz_M = -dz_M; thx_M = -thx_M;
thy_E = My5/lathe.headStockTorsionStiffness;

% Sum contributions due to Forces and Moments
thx_E = thx_F + thx_M;
dz_E = dz_F + dz_M;

% CONSTRUCT HTMS
[HTMs.E, HTMs_d.E, HTMs_m.E] = HTM(lathe.L4, lathe.L5, 0, dx_E, dy_E, dz_E, thx_E, thy_E, thz_E);

%% (F) RAIL SYSTEM
% Assumptions
% - Two point loads act on the rail, representing the rail bearings.
% Average deflection at the two loads is used as the deflection of the
% carriage (x), and the difference over the bearing spacing is used for
% rotation (y).
% - Lead screw is assumed to dominant stiffness in z-direction.

% DIMS AND PROPERTIES (see drawing from Sahil's Notes)
L = lathe.totalRailLength;
a = lathe.distFirstBearing;
b = lathe.railBearingSpacing;
d = lathe.railSpacing;
e = lathe.L8 + lathe.L10;
f = b/2 - lathe.L12;
g = lathe.thirdBearingOffset;
h = lathe.L7 + lathe.L9;

[~,I6] = shapeProps('circle',lathe.d6,[]); % Obtain Inertia of Rail
[A_lScrew,~] = shapeProps('circle',lathe.dlScrew,[]);  % Obtain Area of Leadscrew

% ERRORS DUE TO RAIL BEING BOWED AND AXIS MISALIGNMENT
% x and thy
maxDevX = lathe.railStraightness*L*cos(lathe.railBowOrientation); % Due to bowing
xDev_b1 = maxDevX*(1-4*(a/L - 1/2)^2) + lathe.railMisalignAngleY*a;
xDev_b2 = maxDevX*(1-4*((a + b)/L - 1/2)^2)+ lathe.railMisalignAngleY*(a+b);

dx_assem = xDev_b1 + g/b*(xDev_b2 - xDev_b1);
thy_assem = (xDev_b2 - xDev_b1)/b;

% y and thx
maxDevY = lathe.railStraightness*L*sin(lathe.railBowOrientation);
yDev_b1 = maxDevY*(1-4*(a/L - 1/2)^2)-lathe.railMisalignAngleX*a;
yDev_b2 = maxDevY*(1-4*((a + b)/L - 1/2)^2)-lathe.railMisalignAngleX*(a+b);

dy_assem = yDev_b1 + g/b*(yDev_b2 - yDev_b1);
thx_assem = (yDev_b1 - yDev_b2)/b;

% Vertical Error of Third Bearing (X-Direction not affected due to flexure)
dy_assem3 = lathe.railStraightness*L*sin(lathe.rail2BowOrientation)*...
    (1 - 4*((a+g)/L-1/2)^2) - lathe.rail2MisalignAngleX*(a+g);

% LOADS ON BEARINGS
% x-loads
R1x = Fn*(b - f)/b - Fa*e/b;
R2x = Fn*f/b + Fa*e/b;

% y-loads
K_y = [1 1 1;...
    d d 0;...
    g -(b-g) 0];
F_y = [Ft;...
    Ft*(d-e) + Fn*h;...
    Ft*(g-f) - Fa*h];
R = K_y\F_y;
R1y = R(1); R2y = R(2); R3y = R(3); % Reaction Loads in assumed + y direction

% DEFLECTIONS
% Deflection of Two Bearings on Single Rail
switch lathe.railEndTreatment
    case 'pin-pin'
        % X-Direction Deflection
        [dx_b1, dx_b2] = twoSimplySupported(R1x,R2x,a,b,L,lathe.E6,I6);  

        % Y-Direction Deflection
        [dy_b1, dy_b2] = twoSimplySupported(R1y,R2y,a,b,L,lathe.E6,I6);
        
    case 'fixed-fixed'
        % X-Direction Deflection
        [dx_b1, dx_b2] = twoFixed(R1x,R2x,a,b,L,lathe.E6,I6);  

        % Y-Direction Deflection
        [dy_b1, dy_b2] = twoFixed(R1y,R2y,a,b,L,lathe.E6,I6);
end

% Deflection of Third Bearing
dy_b3 = simplySupported(R3y,a + g,L,lathe.E6,I6);

% Total Deflection and Rotation
dx_F = (dx_b1 + g/b*(dx_b1 + dx_b2)) + dx_assem;
thy_F = (dx_b2-dx_b1)/b + thy_assem;
dy_F = (dy_b1 + g/b*(dy_b1 + dy_b2)) + dy_assem;
thx_F = (dy_b1-dy_b2)/b + thx_assem;
thz_F = (dy_b3 + dy_assem3 -dy_F)/d;

% Z-Stiffness due to axial compression of the lead screw
dz_F = axial(Fa,L,A_lScrew,lathe.E_lScrew);

% CONSTRUCT HTMS
[HTMs.F, HTMs_d.F, HTMs_m.F] = HTM(0, 0, lathe.L6, dx_F, dy_F, dz_F, thx_F, thy_F, thz_F);

%% (G) CARRIAGE
% Assumptions
% 1. The Carraige is treated as rigid. Verify!
dx_G = 0; dy_G = 0; dz_G = 0; thx_G = 0; thy_G = 0; thz_G = 0;

% CONSTRUCT HTMS
[HTMs.G, HTMs_d.G, HTMs_m.G] = HTM(0, lathe.L7, -lathe.L11, dx_G, dy_G, dz_G, thx_G, thy_G, thz_G);

%% (H) CROSSFEED FLEXURE

% Dimensions and Properties
[A8,~] = shapeProps('circle',lathe.XFeedLeadScrewDiameter,[]);  % Obtain Area of Leadscrew

% Calculate Moments at the base of the tool-post
Mz_toolPost = Fn * lathe.L9 - Ft*lathe.L10;
Mx_toolPost = Fa * lathe.L9 - Ft*lathe.L12;
My_toolPost = -Fa * lathe.L10 + Fn*lathe.L12;

% Misalignment Errors
ey = lathe.L8*lathe.xFeedMisalignmentZ;
ez = -lathe.L8*lathe.xFeedMisalignmentY;

% Deflections
dx_H = -Fn/lathe.xFeedStiffnessX; % Due to axial compression
dy_H = -Ft/lathe.flexureYstiffness + ey + lathe.xFeedErrorM_Y;
dz_H = Fa/lathe.flexureZstiffness + ez + lathe.xFeedErrorM_Z;
thx_H = Mx_toolPost/(lathe.flexureYstiffness*(lathe.flexureMountZSpacing)^2) + lathe.xFeedMisalignmentX + lathe.xFeedErrorM_thX;
thy_H = My_toolPost/(lathe.flexureZstiffness*(lathe.flexureMountXSpacing)^2) + lathe.xFeedMisalignmentY + lathe.xFeedErrorM_thY;
thz_H = Mz_toolPost/(lathe.flexureYstiffness*(lathe.flexureMountXSpacing)^2) + lathe.xFeedMisalignmentZ;

[HTMs.H, HTMs_d.H, HTMs_m.H] = HTM(lathe.L8, 0, 0, dx_H, dy_H, dz_H, thx_H, thy_H, thz_H);

%% (I) TOOL POST
% Tool Post is treated as rigid. Verify!
dx_I = 0; dy_I = 0; dz_I = 0; thx_I = 0; thy_I = 0; thz_I = 0;

[HTMs.I, HTMs_d.I, HTMs_m.I] = HTM(0, lathe.L9, -lathe.L12, dx_I, dy_I, dz_I, thx_I, thy_I, thz_I);


%% (J) CUTTING TOOL
% Assume no twisting of the tool

% DIMS
[A10,I10] = shapeProps('rectangle',lathe.w10,lathe.w10); % Obtain Area and Inertia

% DEFLECTIONS
dx_d = -axial(Fn,lathe.L10,A10,lathe.E10); % Due to axial compression
dx_t = thermalExpansion(lathe.al10,lathe.L10,lathe.delT10); % Due to thermal loads
dx_J = dx_t+dx_d; % Total deformation (x-direction)

thx_J = 0;

[dy, thz] = cantilever_F(Ft,lathe.L10,lathe.E10,I10);
dy_J = -dy;
thz_J = -thz;
[dz_J, thy] = cantilever_F(Fa,lathe.L10,lathe.E10,I10);
thy_J = -thy;

% CONSTRUCT HTMS
[HTMs.J, HTMs_d.J, HTMs_m.J] = HTM(lathe.L10,0,0, dx_J, dy_J, dz_J,...
    thx_J, thy_J, thz_J);


%% CALCULATE DEFLECTION OF PART WRT HEADSTOCK
% Deflection of Part with respect to HeadStock
i = 'A':'E';
P_p = [0;0;0;1];
Pd_p = [0;0;0;1];

% Track Contributions from the different elements
nElems = length(i);
elemSensitivity1 = zeros(length(P_p),nElems);
elemSensitivity1(4,:) = 1; 

for j = 1:length(i)
    P_p = HTMs.(i(j))*P_p;
    Pd_p = HTMs_d.(i(j))*Pd_p;
    
    % To determine the contribution from each element, treat it as "stiff"    
    for k = 1:nElems
        if k == j
            elemSensitivity1(:,k) = HTMs.(i(j))*elemSensitivity1(:,k);
        else
            elemSensitivity1(:,k) = HTMs_d.(i(j))*elemSensitivity1(:,k);
        end
    end
end

for k = 1:nElems
    errorConts1(:,k) = (Pd_p(1:3) - P_p(1:3)) - (elemSensitivity1(1:3,k)-P_p(1:3));
end

%% CALCULATE DEFLECTION OF CUTTING TOOL WRT HEADSTOCK
% Deflection of Cutting-Tool with respect to HeadStock
i = 'J':-1:'F';
P_t = [0;0;0;1];
Pd_t = [0;0;0;1];

% Track Contributions from the different elements
nElems = length(i);
elemSensitivity2 = zeros(length(P_t),nElems);
elemSensitivity2(4,:) = 1; 

% Multiply HTM Matrices together
for j = 1:length(i)
    P_t = HTMs.(i(j))*P_t;
    Pd_t = HTMs_d.(i(j))*Pd_t;
    
    % To determine the contribution from an element, treat it as "stiff"    
    for k = 1:nElems
        if k == j
            elemSensitivity2(:,nElems-k+1) = HTMs.(i(j))*elemSensitivity2(:,nElems-k+1);
        else
            elemSensitivity2(:,nElems-k+1) = HTMs_d.(i(j))*elemSensitivity2(:,nElems-k+1);
        end
    end
end

% Calculate the error contributions
for k = 1:nElems
    errorConts2(:,k) = -((Pd_t(1:3) - P_t(1:3)) - (elemSensitivity2(1:3,k)-P_t(1:3)));
end

%% TOTAL DEFLECTION
% Check to see if the strucutral loop is closed
if any(P_p(1:3)-P_t(1:3))>10^-4
    fprintf('\n Warning: Structural loop is not closed.');
end

% Calculate the total deflection
delP = (Pd_p(1:3)-P_p(1:3)) - (Pd_t(1:3)-P_t(1:3));

% Combine Error contributions
errorConts = [errorConts1,errorConts2];
errorContributions = errorConts';

%% EXTRACT DEFLECTION MATRIX
i = 'A':'J';
deflectionMatrix = zeros(length(i),6);

for j = 1:length(i)
    deflectionMatrix(j,:) = [eval(['dx_' i(j)]),...
        eval(['dy_' i(j)]),...
        eval(['dz_' i(j)]),...
        eval(['thx_' i(j)]),...
        eval(['thy_' i(j)]),...
        eval(['thz_' i(j)])];
end

%% DEFORMED AND UNDEFORMED STRUCTURE CALCULATION
% Calculate Deformed Structure from Headstock to Part
i = 'E':-1:'A';
S1 = zeros(4,length(i)+1); % Undeformed Structure
S1(4,:) = 1;
Sd1 = S1; % Deformed structure
for j = 1:length(i)
    S1(:,j+1) = HTMs.(i(j))*S1(:,j);
    Sd1(:,j+1) = HTMs_m.(i(j))*Sd1(:,j);
end

% Calculate Deformed Structure from Headstock to Cutting Tool
i = 'F':'J';
S2 = zeros(4,length(i)+1); % Undeformed Structure
S2(4,:) = 1;
Sd2 = S2; % Deformed structure
for j = 1:length(i)
    S2(:,j+1) = HTMs.(i(j))*S2(:,j);
    Sd2(:,j+1) = HTMs_m.(i(j))*Sd2(:,j);
end


%% ERROR CONTRIBUTION PLOTS
xticks([1 2 3 4 5 6 7 8 9 10]);
    
if plotControl == 1
    close all
    subplot(2,2,1)
    hold on
    bar(errorConts(1,:));
    title('x-Error Contributions');
    ylabel('Error [mm]');
    xticks([1 2 3 4 5 6 7 8 9 10]);
    xticklabels({'Part','Chuck','Spindle','Bearings','Headstock','Rail',...
        'Carriage','Flexure','Toolpost','Tool'});
    xtickangle(45);
    hold off

    subplot(2,2,2)
    hold on
    bar(errorConts(2,:));
    title('y-Error Contributions');
    ylabel('Error [mm]');
    xticks([1 2 3 4 5 6 7 8 9 10]);
    xticklabels({'Part','Chuck','Spindle','Bearings','Headstock','Rail',...
        'Carriage','Flexure','Toolpost','Tool'});
    xtickangle(45);
    hold off

    subplot(2,2,3)
    hold on
    bar(errorConts(3,:));
    title('z-Error Contributions');
    ylabel('Error [mm]');
    xticks([1 2 3 4 5 6 7 8 9 10]);
    xticklabels({'Part','Chuck','Spindle','Bearings','Headstock','Rail',...
        'Carriage','Flexure','Toolpost','Tool'});
    xtickangle(45);
    hold off

    subplot(2,2,4)
    hold on
    plot3(-S1(1,:),S1(3,:),S1(2,:),'LineWidth',2,'Color','k'); % Undeformed Headstock to Part
    plot3(-S2(1,:),S2(3,:),S2(2,:),'LineWidth',2,'Color','k'); % Undeformed Headstock to Cutting Tool
    plot3(-Sd1(1,:),Sd1(3,:),Sd1(2,:),'LineWidth',1,'Color','r'); % Undeformed Headstock to Part
    plot3(-Sd2(1,:),Sd2(3,:),Sd2(2,:),'LineWidth',1,'Color','r'); % Undeformed Headstock to Cutting Tool
    xlabel('x');
    ylabel('z');
    zlabel('y');
    view(135,20);
    hold off
    hold off
end
end

%% SUPPORTING FUNCTIONS
% Area and Second Moment of Inertia Function
function [A,I] = shapeProps(shape,D1,D2)
switch shape
    case 'circle' %D2 can remain null
        A = pi*(D1.^2)/4;
        I = pi()*D1.^4/64;
    case 'rectangle'
        A = D1.*D2;
        I = D1.^3.*D2./12;
    otherwise
        disp('Cross-section is undefined');
end
end

% Deflection due to Axial Loading
function d_s = axial(F, L, A, E)
d_s = (F*L)/(A*E);
end

% Cantilever Deflection due to Point-Force
function [d_s, th_s] = cantilever_F(F, L, E, I)
d_s = (F*L^3)/(3*E*I);
th_s = (F*L^2)/(2*E*I);
end

% Cantilever Deflection due to Point-Moment
function [d_s, th_s] = cantilever_M(M,L,E,I)
d_s = (M*L^2)/(2*E*I);
th_s = (M*L)/(E*I);
end

% Simply-Supported Deflection due to Point Load (Shighley 1015)
function [d_s] = simplySupported(F,a,L,E,I)

d_s = F*(L-a)*a/(6*E*I*L)*(a^2 + (L-a)^2 - L^2);

end

% Deflection due to Two Loads on a Simply Supported Beam 
%(Superposition of Simple Supports - Intermediate Load (Sigley 1015) )
function [ds_1, ds_2] = twoSimplySupported(R1,R2,a,b,L,E,I)
% Dimesions
c = L-a-b;

% Dimensions based on Reference Figure
a1 = a;
b1 = b + c;
a2 = a + b;
b2 = c;

% Deflections
ds_1 = (R1*b1*a1)/(6*E*I*L)*(a1^2 + b1^2 - L^2) + ... % Due to First Load
    (R2*b2*a1)/(6*E*I*L)*(a1^2 + b2^2 - L^2); % Due to Second Load

ds_2 = (R1*a1*(L-a2))/(6*E*I*L)*(a2^2+a1^2-2*L*a2) + ... % Due to First Load
    (R2*b2*a2)/(6*E*I*L)*(a2^2 + b2^2 - L^2); % Due to Second Load

end

% Deflection due to Two Loads on a Cantilevered Beam
%(Superposition of Fixed-Fixed (Sigley 1020) )
function [ds_1, ds_2] = twoFixed(R1,R2,a,b,L,E,I)
% Dimesions
c = L-a-b;

% Dimensions based on Reference Figure
a1 = a;
b1 = b + c;
a2 = a + b;
b2 = c;

% Deflections
ds_1 = (R1*b1^2*a1^2)/(6*E*I*L^3)*(a1*(3*a1+b1)-3*a1*L) + ... % Due to First Load
    (R2*b2^2*a1^2)/(6*E*I*L^3)*(a1*(3*a2+b2)-3*a2*L); % Due to Second Load

ds_2 = (R1*b1^2*(L-a2)^2)/(6*E*I*L^3)*((L-a2)*(3*b1+a1)-3*b1*L) + ... % Due to First Load
    (R2*b2^2*a2^2)/(6*E*I*L^3)*(a2*(3*a2+b2)-3*a2*L); % Due to Second Load
end

% Deflection due to Belt Tension and Cutting Loads on Spindle
function [Rb1, Rb2, d, th] = spindleDeflection(a,b,c,F,M,T,E,I)

% Length
L = a + b + c;

% Solve for reactions on bearings
Rb1 = c/b*F + 1/b*M - T*(a+b)/b;
Rb2 = (c+b)/b*F + 1/b*M - T*a/b;

% Solve for constants by setting up equations to solve
LHS = [1, -1, 0, 0, 0, 0;...
    0, 1, -1, 0, 0, 0;...
    a, -a, 0, 1, -1, 0;...
    0, a+b,-(a+b),0,1,-1;...
    a, 0, 0, 1, 0, 0;...
    0, a+b, 0, 0, 1, 0];

RHS = [-0.5*a^2*Rb1;...
    0.5*(a+b)^2*Rb2;...
    -1/3*a^3*Rb1;...
    1/3*(a+b)^3*Rb2;...
    -1/6*T*a^3;...
    -1/6*T*(a+b)^3-Rb1*(a+b)^2*(-a/3 + b/6)];

consts = LHS\RHS;

% Solve deflections
d = 1/(E*I)*(1/6*T*L^3 + Rb1*(1/6*L^3 - 1/2*a*L^2) - Rb2*(1/6*L^3 - 1/2*(a+b)*L^2) + consts(3)*L + consts(6));
th = 1/(E*I)*(0.5*T*L^2 + Rb1*(1/2*L^2-a*L) - Rb2*(0.5*L^2 - (a+b)*L) + consts(3));

end

% Linear Thermal Expansion
function d_t = thermalExpansion(alpha, L, del_T)
d_t = alpha*L*del_T;
end

% Construct HTMs.
function [H, H_d, H_m] = HTM(Lx,Ly,Lz,dx,dy,dz,thx,thy,thz)
 H = [1, 0, 0, Lx;... % Undeformed HTM
     0, 1, 0, Ly;...
     0, 0, 1, Lz;...
     0, 0, 0, 1];
 
 H_d = [1, -thz, thy, Lx+dx;... % Deformed HTM
       thz, 1, -thx, Ly+dy;...
       thy, thx, 1, Lz+dz;...
     0, 0, 0, 1];
 
 mf = 10; % Magnification of Error
 
 H_m = [1, mf*-thz, mf*thy, Lx+mf*dx;... % Deformed HTM with magnified error
       mf*thz, 1, mf*-thx, Ly+mf*dy;...
       mf*thy, mf*thx, 1, Lz+mf*dz;...
     0, 0, 0, 1];
end
