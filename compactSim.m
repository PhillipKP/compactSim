function [Els, Ennls, Elasso, ab_lassoList, ab_lsList, ...
    F, usedFilterList, abTrue, noiseList] = ...
    compactSim(plotOn, numMeas, mu, noiseInd, filterType, rmsePlotOn)

% Spectral unmixing simulation
%
% This unmixingSim allows us to compare several different strategies for
% selecting the spectral filters.
%
% Phillip K Poon 29 July 2016


% Load the spectral library, ground truth fractional abundance, precomputed
% noise values, and the spectral filters
load('dataForCompactSim.mat','S','specFiltersVis','abReal','noiseMatrix')


% The ground truth fractional abundance for this noise iteration
abTrue = abReal{noiseInd}{1};

% The number of spectral channels and the number of endmembers
numSpecChan = size(S,1);

% The number of endmembers
numSpec = size(S,2);

% The number of spectral filters
numSpecFilt = size(specFiltersVis,2);



% Generate the true mixed spectra
r = S*abTrue';

% Each noiseInd has a seperate list of noise values
noiseList = noiseMatrix(:,noiseInd);

%% Initial variables for speed

ab_lasso = zeros(numSpec, numMeas);
ab_ls = zeros(numSpec, numMeas);
usedFilterList = zeros(1,numMeas);

% Choose a non-adaptive filter

switch filterType
    
    case 'random'
        % Randomly assigns spectral filters
        usedFilterList = randperm(numSpecFilt);
        
    case 'hybrid18'
        % Hybrid PCA-random assignments
        usedFilterList = hybrid18FiltFunc(S, numSpecFilt, specFiltersVis);
        
end

%% Begin measurement loop

for mInd = 1:numMeas
    
    
    % The filter selected for this measurement number
    filterInd = usedFilterList(mInd);
    
    f = specFiltersVis(:, filterInd );
    
    % Update the collection of spectral filters
    F(mInd,:) = f';
    
    
    if plotOn == 1
        figure(450);
        imagesc(F)
        title('F')
        ylabel('Measurement Number','fontsize',14)
        xlabel('Spectral Channel','fontsize',14)
        set(gca,'FontSize',14)
        
        
    end
    
    
    %% Update H
    H = F*S;
    
    
    %% Actual simulated measurements with noise
    
    g(mInd,:) = f.' * r +  noiseList(mInd);
    
    
    %% Compute the least squares estimate
    
    ab_ls = pinv(H'*H)*H'*g;
    
    % Compute the RMSE of the LS estimate
    
    Els(mInd) = mse_func(ab_ls,abTrue');
    
    
    % Store the estimated fractional abundance
    ab_lsList(:,mInd) = ab_ls;
    
    
    %% Compute the nonnegative LS (NNLS) estimate of the abundance
    % fraction for this pixel
    
    ab_nnls = lsqnonneg(H, g);
    
    
    % Compute the RMSE of the LS estimate
    Ennls(mInd) = mse_func(ab_nnls,abTrue');
    
    %% The unconstrained LASSO for this pixel
    
    if mInd > 1
        
        lambda = mu;
        
        ab_lasso = lasso(H, g, 'lambda', lambda);
        
    elseif mInd == 1
        
        ab_lasso = ab_nnls;
        
    end
    
    % Enforce Positivity
    ab_lasso = ab_lasso.*(ab_lasso > 0);
    
    % Store the estimated fractional abundance
    ab_lassoList(:,mInd) = ab_lasso;
    
    %% Compute metrics
    
    % Compute the RMSE of the LASSO estimate
    Elasso(mInd) = mse_func(ab_lasso, abTrue');
    
    
    
    %% Plots the fractional abundances vs measurements
    
    if plotOn == 1
        
        figure(400);
        subplot(3,1,1);
        plot(ab_lassoList','linewidth',3);
        ylim([0 1]);
        title('Randomly Selected. LASSO');
        
        subplot(3,1,2);
        plot(ab_lsList','linewidth',3);
        ylim([0 1]);
        title('Randomly Selected. LS');
        
        subplot(3,1,3);
        plot(repmat(abTrue,[40 1]),'linewidth',3);
        ylim([0 1]);
        title('Truth');
        drawnow
        
    end
    
    if rmsePlotOn == 1 && mInd > 2
        
        figure(402);
        semilogy(Elasso);
        xlim([2 mInd])
        legend('LASSO RMSE')
        drawnow
        
    end
    
    
end


end
