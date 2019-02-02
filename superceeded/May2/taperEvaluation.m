clear all;
clc;

myLathe = lathe;

%% Initialization
errorMatrix = struct;

% Part Details
L2Dratio = 4;
maxLength = 80;
partDia = 2:4:30; %[mm]

% Base Loads
Fn = 12;%[N] Steel: 12, Aluminum: 5;
Ft = 15;%[N] Steel: 15, Aluminum: 6;
Fa = 12;%[N] Steel: 12, Aluminum: 5;
loadVarFactor = 0.50; % +/- on base load

%% CALCULATE ERROR

for i = 1:length(partDia)
    myLathe.D1 = partDia(i);
    errorMatrix.case(i).partDiameter = partDia(i);
    
    partLengths = 1:((partDia(i)*L2Dratio)/20):min((partDia(i)*L2Dratio),maxLength);
    errorMatrix.case(i).partLengths = partLengths';
    
    for j = 1:length(partLengths)
        
        myLathe.L1 = partLengths(j);
        
        % At Base Load
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
        errorMatrix.case(i).baseLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).baseLoad.individualConts(j,:) = errorConts(:,1)';
        
        % At Lower Load
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1-loadVarFactor)*Fn, (1-loadVarFactor)*Ft, (1-loadVarFactor)*Fa, myLathe, 0);
        errorMatrix.case(i).lowerLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).lowerLoad.individualConts(j,:) = errorConts(:,1)';
        
        % At Upper Load
        [delP, errorConts, ~, ~, ~] = HTM_Error_Estimator((1+loadVarFactor)*Fn, (1+loadVarFactor)*Ft, (1+loadVarFactor)*Fa, myLathe, 0);
        errorMatrix.case(i).upperLoad.totalErrors(j,:) = delP';
        errorMatrix.case(i).upperLoad.individualConts(j,:) = errorConts(:,1)';
        
        fprintf('\nRun %3.0f of %4.0f\n',(i-1)*length(partLengths)+j,length(partDia)*length(partLengths));
    end
end

%% COMPARE TOTAL DIAMETER ERROR OVER DIFFERENT PART LENGTHS
hold on
for i = 1:length(errorMatrix.case)
    plot(errorMatrix.case(i).partLengths(:,1),2*1000*errorMatrix.case(i).baseLoad.totalErrors(:,1));
    legendInfo{i} = num2str(errorMatrix.case(i).partDiameter);
end
ax = gca;
xticks('auto');
xlabel('Part Length [mm]');
ylabel('Modeled Diametrical Error [um]');
xlim([0,100]);
legend(legendInfo);
hold off


%% COMPARE ERROR SOURCES FOR A FIXED PART DIAMETER
figure()
caseNum = 2;
hold on
for i = 1:length(errorMatrix.case(caseNum).baseLoad.individualConts(1,:))
    plot(errorMatrix.case(caseNum).partLengths(:,1),2*1000*errorMatrix.case(caseNum).baseLoad.individualConts(:,i));
end
legend('Part','Chuck','Spindle','Bearings','Headstock','Rail',...
        'Carriage','Flexure','Toolpost','Tool');
ax = gca;
xticks('auto');
xlabel('Part Length [mm]');
ylabel('Modeled Contribution to Diameter Error [um]');
hold off

%% POST PROCESS ERROR CHART
compContributions = zeros(length(errorMatrix.case),1 + 10);

for i = 1:length(errorMatrix.case)
    compContributions(i,1) = errorMatrix.case(i).partDiameter;
    compContributions(i,2:end) = 2*(max(errorMatrix.case(i).upperLoad.individualConts)-...
        min(errorMatrix.case(i).lowerLoad.individualConts))*1000;
end

%% OUTPUT TO EXCEL FILE
filename = 'errorBudget.xlsx';
headers = {'Diameter [mm]','Part','Chuck','Spindle','Bearings','Headstock','Rail',...
         'Carriage','Flexure','Toolpost','Tool'};
xlswrite(filename,compContributions,'TaperAnalysis','B3');
xlswrite(filename,headers,'TaperAnalysis','B2');
