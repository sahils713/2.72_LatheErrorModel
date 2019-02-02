clear all;
clc;

myLathe = lathe;
myLatheCopy = lathe;

%% Initialization
errorMatrix = struct;

% Part Details
L2Dratio = 4;
maxDiameter = 30; %[mm]
partLength = 5:5:80; %[mm]

% Base Loads
Fn = 12;%[N] Steel: 12, Aluminum: 5;
Ft = 15;%[N] Steel: 15, Aluminum: 6;
Fa = 12;%[N] Steel: 12, Aluminum: 5;
loadVarFactor = 0.50; % +/- on base load

%% CALCULATE ERROR

for i = 1:length(partLength)
    myLathe.L1 = partLength(i);
    errorMatrix.case(i).partLength = partLength(i);
    
    partDiameters = (partLength(i)/L2Dratio):1:maxDiameter;
    errorMatrix.case(i).partDiameters = partDiameters';
    
    for j = 1:length(partDiameters)
        
        myLathe.D1 = partDiameters(j);
        
        % At Base Load with Error Motions
        myLathe.xFeedErrorM_thX = myLatheCopy.xFeedErrorM_thX;
        myLathe.xFeedErrorM_thY = myLatheCopy.xFeedErrorM_thY;
        myLathe.xFeedErrorM_Z = myLatheCopy.xFeedErrorM_Z;
        myLathe.xFeedErrorM_Y = myLatheCopy.xFeedErrorM_Y;
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
        errorMatrix.case(i).baseLoad_EM.totalErrors(j,:) = delP';
        errorMatrix.case(i).baseLoad_EM.individualConts(j,:) = errorConts(:,3)';
        
        % At Base Load without Error Motions
        % Set all error motions to zero
        myLathe.xFeedErrorM_thX = 0;
        myLathe.xFeedErrorM_thY = 0;
        myLathe.xFeedErrorM_Z = 0;
        myLathe.xFeedErrorM_Y = 0;
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
        errorMatrix.case(i).baseLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).baseLoad.individualConts(j,:) = errorConts(:,3)';
        
        % At Lower Load (without Error Motions)
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1-loadVarFactor)*Fn, (1-loadVarFactor)*Ft, (1-loadVarFactor)*Fa, myLathe, 0);
        errorMatrix.case(i).lowerLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).lowerLoad.individualConts(j,:) = errorConts(:,3)';
        
        % At Upper Load (without Error Motions)
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1+loadVarFactor)*Fn, (1+loadVarFactor)*Ft, (1+loadVarFactor)*Fa, myLathe, 0);
        errorMatrix.case(i).upperLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).upperLoad.individualConts(j,:) = errorConts(:,3)';
        
        fprintf('\nRun %3.0f of %4.0f\n',(i-1)*length(partLength) + j,length(partDiameters)*length(partLength));
    end
end

%% COMPARE TOTAL DIAMETER ERROR OVER DIFFERENT PART LENGTHS
hold on
for i = 1:length(errorMatrix.case)
    plot(errorMatrix.case(i).partDiameters(:,1),errorMatrix.case(i).baseLoad_EM.totalErrors(:,3));
    legendInfo{i} = num2str(errorMatrix.case(i).partLength);
end
ax = gca;
xticks('auto');
xlabel('Part Diameter [mm]');
ylabel('Modeled Length Error [um]');
legend(legendInfo);
hold off


%% COMPARE ERROR SOURCES FOR A FIXED PART DIAMETER
figure()
caseNum = 2;
hold on
for i = 1:length(errorMatrix.case(caseNum).baseLoad_EM.individualConts(1,:))
    plot(errorMatrix.case(caseNum).partDiameters(:,1),1000*errorMatrix.case(caseNum).baseLoad_EM.individualConts(:,i));
end
legend('Part','Chuck','Spindle','Bearings','Headstock','Rail',...
        'Carriage','Flexure','Toolpost','Tool');
ax = gca;
xticks('auto');
xlabel('Part Diameter [mm]');
ylabel('Modeled Contribution to Length Error [um]');
hold off

%% POST PROCESS ERROR CHART
compContributions = zeros(length(errorMatrix.case),1 + 10);

for i = 1:length(errorMatrix.case)
    compContributions(i,1) = errorMatrix.case(i).partLength;
    % Add the error contribution from error motions and load variations (stiffness)
    compContributions(i,2:end) =...
        (max(errorMatrix.case(i).baseLoad_EM.individualConts - errorMatrix.case(i).baseLoad.individualConts) +...
        (max(errorMatrix.case(i).upperLoad.individualConts)-...
        min(errorMatrix.case(i).lowerLoad.individualConts)))*1000;
end

%% OUTPUT TO EXCEL FILE
filename = 'errorBudget.xlsx';
headers = {'Length [mm]','Part','Chuck','Spindle','Bearings','Headstock','Rail',...
         'Carriage','Flexure','Toolpost','Tool'};
xlswrite(filename,compContributions,'FlatnessAnalysis','B3');
xlswrite(filename,headers,'FlatnessAnalysis','B2');
