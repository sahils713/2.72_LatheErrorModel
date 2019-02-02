clear all;
clc;

myLathe = lathe;
myLatheCopy = lathe;

%% Initialization
errorMatrix = struct;

% Part Details
L2Dratio = 4;
partDia = 2:2:30; %[mm]
partLength = 4:4:80; %[mm]
flatnessError = zeros(length(partLength),length(partDia));

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

%% CALCULATE ERROR

for i = 1:length(partLength)
    myLathe.L1 = partLength(i);
    errorMatrix.case(i).partLength = partLength(i);
      
    for j = 1:length(partDia)
        
        if (partLength(i)/partDia(j))>L2Dratio
            errorMatrix.case(i).totalErrors(j,:) = [NaN;NaN;NaN];
            flatnessError(i,j) = NaN;
        else
            errorMatrix.case(i).partDiameters(j) = partDia(j);
            myLathe.D1 = partDia(j);
            cuttingDepths = linspace(1,partDia(j)/2,5);

            for k = 1:length(cuttingDepths)

                myLathe.cuttingDepth = cuttingDepths(k);

                % At Base Load with Error Motions
                myLathe.xFeedErrorM_thX = myLatheCopy.xFeedErrorM_thX;
                myLathe.xFeedErrorM_thY = myLatheCopy.xFeedErrorM_thY;
                myLathe.xFeedErrorM_Z = myLatheCopy.xFeedErrorM_Z;
                myLathe.xFeedErrorM_Y = myLatheCopy.xFeedErrorM_Y;
                [delP, ~, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
                baseLoad_EMErrors(k,:) = delP';
                %errorMatrix.case(i).baseLoad_EM.individualConts(j,:) = errorConts(:,3)';

                % At Base Load without Error Motions
                % Set all error motions to zero
                myLathe.xFeedErrorM_thX = 0;
                myLathe.xFeedErrorM_thY = 0;
                myLathe.xFeedErrorM_Z = 0;
                myLathe.xFeedErrorM_Y = 0;
                [delP, ~, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
                baseLoadErrors(k,:) = delP';

                % At Lower Load (without Error Motions)
                [delP, ~, ~, ~, ~] = HTM_Error_Estimator((1-loadVarFactor)*Fn, (1-loadVarFactor)*Ft, (1-loadVarFactor)*Fa, myLathe, 0);
                lowerLoadErrors(k,:) = delP';

                % At Upper Load (without Error Motions)
                [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1+loadVarFactor)*Fn, (1+loadVarFactor)*Ft, (1+loadVarFactor)*Fa, myLathe, 0);
                upperLoadErrors(k,:) = delP';
            end
            errorMatrix.case(i).totalErrors(j,:) = ...
                (max(upperLoadErrors - lowerLoadErrors) +...
            (max(baseLoad_EMErrors - baseLoadErrors)) + ...
            (max(baseLoadErrors)-min(baseLoadErrors)))*1000;

            flatnessError(i,j) = errorMatrix.case(i).totalErrors(j,3);
        end
        
        fprintf('\nRun %3.0f of %4.0f\n',(i-1)*length(partLength) + j,length(partDia)*length(partLength));
    end
end


%% OUTPUT TO EXCEL FILE
filename = 'errorBudget.xlsx';
flatnessErrorTable = [0,partDia;partLength',flatnessError];
xlswrite(filename,flatnessErrorTable,'FlatnessAnalysis','C4');