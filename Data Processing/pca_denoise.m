function [denoised,nPCs,eig,expVar,mu] = pca_denoise(data,nPCs)
% This function removes noise by removing a given number of
% Principal Components before reconstructing the data.
%
% Renan M. Costa
%
% July 5th 2020, updated to automatically find optimal number of PCs to
% remove if a single input argument is provided. If two input arguments are
% provided, the second argument is the number of PCs to remove. Output data
% ("denoised") is now z-scored.

%% Parameters

plotZCovs = 0; % Whether to plot covariance of z-scored data as a function of PCs removed.

%% Run PCA
dataZsc = zscore(data); %Z-score input data.
[eig, Z,latent,tsquared,expVar,mu] = pca(dataZsc);

% The projection onto principal component space ("Z") can be obtained by
% multiplying the mean-subtracted data by eig. That is:
% Z = (data - repmat(mu,size(data,1),1))*eig;

%% Denoise
if nargin~=1 % If a number of PCs is given as an input, remove that many PCs.
    eigZeroed = eig;
    eigZeroed(:,1:nPCs) = 0;
    
    % Reconstruct centered data without removed components.
    denoised = Z * eigZeroed';
    
    % Complete reconstruction by adding back the means of each trace, which the
    % PCA algorithm subtracts when centering the data.
    denoised = denoised + repmat(mu,size(Z,1),1);
    denoised = zscore(denoised);
    
else % If no number of PCs is given, find the optimal number to remove. 
     % This is done by minimizing the total covariance of the z-scored 
     % denoised data. Note that the covariance matrix of the z-scored data 
     % is equivalent to the correlation matrix of the data.
    
    zCovs = nan(size(dataZsc,2),1);
    zCovs(1) = mean(abs(cov(dataZsc)),'all'); % Mean absolute covariance of z-scored data.
    for A = 1:size(dataZsc,2)-1 % Iterate through removing all numbers of PCs, and determine resulting total covariance in each case.
        eigZeroed = eig;
        eigZeroed(:,1:A) = 0;
        
        % Reconstruct centered data without removed components.
        denoised = Z * eigZeroed';
        
        % Complete reconstruction by adding back the means of each trace, that the
        % PCA algorithm subtracts when centering the data.
        denoised = denoised + repmat(mu,size(Z,1),1);
        denoised = zscore(denoised);
        
        zCovs(A+1) = mean(abs(cov(denoised)),'all'); % Mean absolute covariance of z-scored denoised data.
    end
    % Find minimum covariance and number of PCs to remove.
    [~,idx] = min(zCovs);
    nPCs = idx - 1; % Need to subtract one because the first element of zCovs corresponds to removal of 0 PCs.
    
    % Recursive step. Call function again with number of PCs above
    % as an input.
    [denoised,nPCs,eig,expVar,mu] = pca_denoise(data,nPCs);
end

%% Plot covariance of z-scored data as a function of PCs removed.

if plotZCovs && exist('zCovs','var')
    f1 = figure;
    ax1 = axes; hold on;
    
    hline = plot(0:size(zCovs,1)-1,zCovs,'LineWidth',2);
    hsc = scatter(nPCs,min(zCovs),[],'k','filled');
    
    %ax.YLim(1) = 0; % Commented out, only useful in linear scale.
    set(ax1,'YScale','log')
    legend(hsc,'Minimum','Location','best')
    ylabel 'Mean Absolute Correlation'
    xlabel 'Principal Components Removed'
    title 'Correlation'
    ax1.TickDir = 'out';
    
end

end