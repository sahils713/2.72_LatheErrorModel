myLathe = lathe;

%% Initialization
cuttingToolOffset = 0:1:40; %[mm]
errorX = zeros(length(cuttingToolOffset),3); % Error Terms

Fn = 5;
Ft = 6;
Fa = 5;

%% CALCULATE ERROR

for i = 1:length(cuttingToolOffset)
    myLathe.L12 = cuttingToolOffset(i);
    [delP, errorCont, ~, ~, ~] = HTM_Error_Estimator(Fn, Ft, Fa, myLathe, 0);
    
    errorX(i,:) = errorCont(8,:);
        
    fprintf('\nRun %3.0f of %4.0f\n',i,length(cuttingToolOffset));
end

%% PLOT ERROR
hold on
plot(cuttingToolOffset,errorX(:,1)*1000);
plot(cuttingToolOffset,errorX(:,2)*1000);
plot(cuttingToolOffset,errorX(:,3)*1000);
legend('x','y','z');
ax = gca;
xticks('auto');
xlabel('Cutting Tool Offset [mm]');
ylabel('Modeled Error [um]');
set(gca,'FontSize',14);
hold off
