partDia = 10:1:13; %[mm]
partLength = 20:5:100; %[mm]

errorX = size(length(partDia),length(partLength)); % Error Terms

Fn = 10;
Ft = 100;
Fa = 10;

%% CALCULATE ERROR

for i = 1:length(partDia)
    lathe.D1 = partDia(i);
    
    for j = 1:length(partLength);
        
        lathe.L1 = partLength(j);
        [delP, ~, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, lathe, 0);
        
        errorX(i,j) = delP(2);
        
    end
    fprintf('\nRun %3.0f of %4.0f',i*j,length(partDia)*length(partLength));
end

%% PLOT ERROR
figure()
hold on
surf(partLength,partDia,errorX.*1000);
xlabel('Part Length [mm]');
ylabel('Part Diameter [mm]');
zlabel('Modeled Error [um]');
hold off
