clear;
clc;

myLathe = lathe;

%% INITIALIZATION

% Part Details
L2Dratio = 4; % Maximum Length to Width Ratio
partDia = 2:2:30; %[mm]
partLength = 4:4:80; %[mm]

% Error Matrix
lowLoadError = zeros(length(partDia),length(partLength),3);
baseLoadError = zeros(length(partDia),length(partLength),3);
highLoadError = zeros(length(partDia),length(partLength),3);

% Base Loads
if myLathe.mat1 == 'St'
    Fn = 12;%[N] Steel: 12, Aluminum: 5;
    Ft = 15;%[N] Steel: 15, Aluminum: 6;
    Fa = 12;%[N] Steel: 12, Aluminum: 5;
elseif myLathe.mat1 == 'Al'
    Fn = 5;%[N] Steel: 12, Aluminum: 5;
    Ft = 6;%[N] Steel: 15, Aluminum: 6;
    Fa = 5;%[N] Steel: 12, Aluminum: 5;
else
    disp('Part material specification error');
end
loadVarFactor = 0.50; % +/- on base load

% Base temperatures
delT1 = 8;
delT2 = 8;
delT3 = 10;
delT10 = 10;
tempVar = 3; % +/- deg C on base temperature

%% CALCULATE ERRORS

for i = 1:length(partDia)
    myLathe.D1 = partDia(i);
    
    for j = 1:length(partLength)
        
        if (partLength(j)/partDia(i))>L2Dratio
            lowLoadError(i,j,1:3) = NaN;
            baseLoadError(i,j,1:3) = NaN;
            highLoadError(i,j,1:3) = NaN;
        else
            myLathe.L1 = partLength(j);
            
            %At Base Load
            myLathe.delT1 = delT1; myLathe.delT2 = delT2;
            myLathe.delT3 = delT3; myLathe.delT10 = delT10;
            [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
            baseLoadError(i,j,1:3) = delP;
            errorContsBaseLoad(i,j,:) = errorConts(:,1);
            
            
            % At Lower Load and Lower Temperatures
            myLathe.delT1 = delT1-tempVar; myLathe.delT2 = delT2-tempVar;
            myLathe.delT3 = delT3-tempVar; myLathe.delT10 = delT10-tempVar;
            [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1-loadVarFactor)*Fn, (1-loadVarFactor)*Ft, (1-loadVarFactor)*Fa, myLathe, 0);
            lowLoadError(i,j,1:3) = delP;
            errorContsLowLoad(i,j,:) = errorConts(:,1);

            % At Upper Load and Higher Temperatures
            myLathe.delT1 = delT1+tempVar; myLathe.delT2 = delT2+tempVar;
            myLathe.delT3 = delT3+tempVar; myLathe.delT10 = delT10+tempVar;
            [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1+loadVarFactor)*Fn, (1+loadVarFactor)*Ft, (1+loadVarFactor)*Fa, myLathe, 0);
            highLoadError(i,j,1:3) = delP;
            errorContsHighLoad(i,j,:) = errorConts(:,1);
        end
              
        fprintf('\nRun %3.0f of %4.0f\n',(i-1)*length(partLength)+j,length(partDia)*length(partLength));
    end
end

%% RECHUCKING ERROR
% Rechucking error was measured to be 0.003" on the diameter on a gauge 
% pin 2" away from the tip of the chuck. If parts are being turned down
% more than 0.003", this error should not affect the repeatability.

%% POST-PROCESS REPEATABILITY ERROR
% Multiplied by 2 to change offset error into a diametrical error
% Mutiplied by 1000 to change from mm to um
repeatabilityError = 2*1000*(highLoadError(:,:,1)-lowLoadError(:,:,1));
repeatabilityErrorConts = 2*1000*(errorContsHighLoad - errorContsLowLoad);

%% PLOT ERROR AND REPEATABILITY

% Base Load Error
hold on
surf(partLength,partDia,baseLoadError(:,:,1)*1000);
surf(partLength,partDia,lowLoadError(:,:,1)*1000);
surf(partLength,partDia,highLoadError(:,:,1)*1000);
xticks('auto');
xlabel('Part Length [mm]');
ylabel('Part Diameter [mm]');
zlabel('Modeled x-Error [um]');
hold off

figure()
% Repeatability
hold on
surf(partLength,partDia,repeatabilityError);
xticks('auto');
xlabel('Part Length [mm]');
ylabel('Part Diameter [mm]');
zlabel('Repeatability Error [um]');
hold off

%% OUTPUT TO EXCEL FILE
filename = 'errorBudget.xlsx';
repeatabilityErrorTable = [0,partLength;partDia',repeatabilityError];
xlswrite(filename,repeatabilityErrorTable,'RepeatabilityAnalysis','C4');
