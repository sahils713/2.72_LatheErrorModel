myLathe = lathe;

%% Initialization
partDia = 0:1:32; %[mm]
partLength = 0:2:100; %[mm]

errorX = size(length(partDia),length(partLength)); % Error Terms

Fn = 5;
Ft = 6;
Fa = 5;

%% CALCULATE ERROR

for i = 1:length(partDia)
    myLathe.D1 = partDia(i);
    
    for j = 1:length(partLength)
        
        myLathe.L1 = partLength(j);
        [~, errorCont, ~, ~, ~] = HTM_Error_Estimator(1.5*Fn, 1.5*Ft, 1.5*Fa, myLathe, 0);
        
        errorX(i,j) = 2*errorCont(1,1); % Convert to diametrical error
        
        fprintf('\nRun %3.0f of %4.0f\n',(i-1)*length(partLength)+j,length(partDia)*length(partLength));
      
    end
end

%% PLOT ERROR
hold on
surf(partLength,partDia,errorX.*1000);
ax = gca;
xticks('auto');
xlabel('Part Length [mm]');
ylabel('Part Diameter [mm]');
zlabel('Modeled Error [um]');
zlim([0,25]);
caxis([0,25]);
h = colorbar;
ylabel(h,'Modeled Error [um]');
set(h,'FontSize',14);
set(gca,'FontSize',14);
hold off
