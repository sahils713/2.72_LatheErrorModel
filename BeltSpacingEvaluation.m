myLathe = lathe;
myLathe.D1 = 20;

%% Initialization
distanceBelt = 10:1:40; %[mm]
partLength = [80]; %[mm]

errorX = size(length(distanceBelt)); % Error Terms

Fn = 10;
Ft = 20;
Fa = 10;

%% CALCULATE ERROR
for j = 1:length(partLength)
    myLathe.L1 = partLength(j);
    
    for i = 1:length(distanceBelt)
    myLathe.beltDistance = distanceBelt(i);
    
    [delP, ~, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
    errorY(i) = delP(2);
   fprintf('\nRun %3.0f of %4.0f\n',i*j,length(distanceBelt)*length(partLength));
    end
    
    %% PLOT ERROR
    hold on
    plot(distanceBelt,errorY.*1000);

end
ax = gca;
xticks('auto');
xlabel('Belt Spacing from Nearest Bearing [mm]');
ylabel('Error [um]');
