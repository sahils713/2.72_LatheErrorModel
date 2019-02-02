myLathe = lathe;
myLathe.D1 = 20;

%% Initialization
distanceBearings = 10:1:160; %[mm]
partLength = [80]; %[mm]

errorX = size(length(distanceBearings)); % Error Terms

Fn = 5;
Ft = 6;
Fa = 5;

%% CALCULATE ERROR
for j = 1:length(partLength)
    myLathe.L1 = partLength(j);
    
    for i = 1:length(distanceBearings)
    myLathe.shaftBearingSpacing = distanceBearings(i);
    
    [delP, ~, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
    errorX(i) = delP(1);
   fprintf('\nRun %3.0f of %4.0f\n',i*j,length(distanceBearings)*length(partLength));
    end
    
    %% PLOT ERROR
    hold on
    plot(distanceBearings,errorX.*1000);

end
ax = gca;
xticks('auto');
xlabel('Bearing Spacing [mm]');
ylabel('Error [um]');
