function [] = launchSim()
% [] = launchSim()
%
% This version of the unmixing simulation varies mu (aka lambda) instead of
% c. Before when we studied how the unmixing performance changes with c, c
% was used to compute mu (aka lambda) but that meant that mu would change
% with every measurement.
% So in this version of the simulation we want to vary mu from 0.1 to 0.7
% in 100 steps and see if we see any different.

%% User Defined Variables

% Number of measurements
numMeas = 20;

% Signal to Noise
snr = 2;

adaptiveOn = 0;

muList = [0.006];

noiseIterList = 1

plotOn = 0

filterType = 'hybrid18'

rmsePlotOn = 1
%%


for mu = muList
    
    ElsList = [];
    EnnlsList = [];
    ElassoList = [];
    ElassoMinList = [];
    
    for noiseInd = noiseIterList
        
        
        [Els, Ennls, Elasso, ab_lassoList, ab_lsList, ...
            F, usedFilterList, abTrue, noiseList] = ...
            compactSim(plotOn, numMeas, mu, noiseInd, filterType, rmsePlotOn);
        
        ElsList = [ElsList; Els];
        EnnlsList = [EnnlsList; Ennls];
        ElassoList = [ElassoList; Elasso];
        ElassoMinList = [ElassoMinList min(Elasso)]
        
        
        
        ElsAvg = mean(ElsList,1);
        EnnlsAvg = mean(EnnlsList,1);
        ElassoAvg = mean(ElassoList,1);
        
        % What the saved file will be called
        fileName = ['result_noiseInd_' num2str(noiseInd) '.mat'];
        
        % Actually saves the file
        
        save(fileName,...
            'ElsAvg','ElassoAvg','EnnlsAvg', 'snr','mu','numMeas',...
            'noiseIterList','muList','snr','noiseInd');
        
        
    end
    
    
    figure(500);
    semilogy(ElsAvg,'linewidth',3)
    hold all
    semilogy(ElassoAvg, 'linewidth',3)
    legend('LS','LASSO')
    xlabel('Measurements')
    ylabel('RMSE')
    title(['SNR = ' num2str(10^snr) ', mu = ' num2str(mu) ''])
    set(gca,'FontSize',14)
    grid on
    ylim([10^-4 10^2])
    hold off
    drawnow
    

    
    
end


end
