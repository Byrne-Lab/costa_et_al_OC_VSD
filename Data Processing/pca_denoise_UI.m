function pca_denoise_UI
%% Parameters

% Load data. .
load("E:\Renan\Operant Conditioning\Test Data\SpikeAmpTestData_16-Jul-2020.mat",'vsd_filtered_ns','vsddata');
%[b,a] = ellip(2,0.1,40,[10 100]*2/1000);
%vsd_ellip = filtfilt(b,a,vsddata);
data = zscore(vsd_filtered_ns);
%data = zscore(vsd_ellip);
%   VARIABLES                   DESCRIPTION
%   data                        Matrix containing VSD traces. Rows are timepoints and columns are neurons.

% Acquisition rate of data in Hz;
acqRate = 1000;

% Whether or not to use variance threshold. 1 uses variance threshold, 0
% finds optimal number of PCs to remove.
useThr = 0;

% Variance threshold. Principal components adding up to this variance explained will
% be removed.
varThr = 85;

%% Run PCA
[eigV, Z,~,~,expVar,mu] = pca(data);
%data = zscore(data);
%   OUTPUT VARIABLES            DESCRIPTION
%   eigV                        "Components", or eigenvectors. Frequently called weights or loadings.
%   Z                           "Projection", or projection of data onto principal component space. Frequently called scores.
%   expVar                      Amount of variance explained by each principal component.
%   mu                          Mean of each column of "data". The PCA algorithm subtracts the mean of each trace, so the mean needs to be added back to reconstruct the traces later.

%% Denoise traces

% Remove PCs
if useThr
    nPCs = find(cumsum(expVar)>varThr,1);
    [denoised,~] = removeNoise(nPCs);
else
    [denoised,nPCs] = removeNoise();
end

varThr = sum(expVar(1:nPCs));
% 
% % Zero out components up to the number above
% eigZeroed = eigV;
% eigZeroed(:,1:nPCs) = 0;
% 
% % Reconstruct centered data without removed components.
% denoised = Z * eigZeroed';
% 
% % Complete reconstruction by adding back the means of each trace, that the
% % PCA algorithm subtracts when centering the data.
% denoised = denoised + repmat(mu,size(Z,1),1);

% Note:
% The projection onto principal component space ("Z") can be obtained by
% multiplying the mean-subtracted data by eig. That is:
% Z = (data - repmat(mu,size(data,1),1))*eig;

%% Detect spikes

[orig_spiketimes,orig_spTraces,~,~]=detect_spikes_vsd(data);
[denoised_spiketimes,denoised_spTraces,~,~]=detect_spikes_vsd(denoised);

% %% Reorder by firing rate
%
% [~,idx] = sort(sum(denoised_spiketimes~=0,1),'descend');
%
% data = data(:,idx);
% denoised = denoised(:,idx);
% orig_spiketimes = orig_spiketimes(:,idx);
% denoised_spiketimes = denoised_spiketimes(:,idx);
% eigV = eigV(idx,:);
% mu = mu(idx);

%% Obtain individual and average spike traces

% Compute average spike traces
orig_avgSp = cellfun(@(x) mean(x,1),orig_spTraces,'uni',0);
orig_avgSp = cell2mat(orig_avgSp')'; % Double transpose necessary for cell2mat.
denoised_avgSp = cellfun(@(x) mean(x,1),denoised_spTraces,'uni',0);
denoised_avgSp = cell2mat(denoised_avgSp')';

%% Plotting

% Parameters

neurPlot = 1:5:50; % Specify which neurons to plot.

%% Initiate figure
if ~isempty(findobj('Type','Figure','Name','Denoising Adjustment'))
    close('Denoising Adjustment')
end
f1 = figure('Name','Denoising Adjustment');
f1.Position = [340   50   1200   920];

%% Plot neuron traces

% Offset to add to spiketimes for plotting.
% (Spiketimes are at onset of each spike. This adds a few
% datapoints to approximate the peak.)
spTimeOffset = 0;

% Preallocate all graphic objects
[ax1, ax2, neurLabel,...
    hOrig, hDenoised, ...
    hScOrig, hScDenoised, ...
    hAvgOrig, hAvgDenoised] = ...
    deal(gobjects(length(neurPlot),1));

for A = 1:length(neurPlot)
    
    % Crate axes for main traces
    ax1(A) = axes; hold on;
    ax1(A).Position = [0.025 0.95-0.4/length(neurPlot)*A 0.90 0.4/length(neurPlot)];
    
    % Label each neuron
    neurLabel(A) = annotation('textbox',...
        [0.0045 0.95-0.4/length(neurPlot)*A 0.019 0.4/length(neurPlot)],....
        'String',neurPlot(A),'EdgeColor','none','HorizontalAlignment','right',...
        'VerticalAlignment','middle','FontSize',12,'Margin',0);
    
    % Plot main traces
    hOrig(A) = plot(data(:,neurPlot(A)));
    hDenoised(A) = plot(denoised(:,neurPlot(A)));
    
    % Prepare spike times
    tempSTO = orig_spiketimes(:,neurPlot(A));
    tempSTO(tempSTO==0) = [];
    tempSTO = tempSTO + spTimeOffset;
    tempSTD = denoised_spiketimes(:,neurPlot(A));
    tempSTD(tempSTD==0) = [];
    tempSTD = tempSTD + spTimeOffset;
    
    % Plot spike times
    try
        hScOrig(A) = plot(tempSTO,data(tempSTO,neurPlot(A)),'|','LineWidth',2);
    catch
        hScOrig(A) = plot(NaN,NaN,'|'); % If plotting fails (e.g., there are no spikes), still create a non-empty line object to be modified by nested functions.
    end
    
    try
        hScDenoised(A) = plot(tempSTD,denoised(tempSTD,neurPlot(A)),'|','LineWidth',2);
    catch
        hScDenoised(A) = plot(NaN,NaN,'|');
    end
    
    % Create axes for average spike traces
    ax2(A) = axes; hold on;
    ax2(A).Position = [0.94 0.95-0.4/length(neurPlot)*A 0.05 0.4/length(neurPlot)];
    
    % Plot average traces
    hAvgOrig(A) = plot(orig_avgSp(:,neurPlot(A)),'LineWidth',1);
    hAvgDenoised(A) = plot(denoised_avgSp(:,neurPlot(A)),'LineWidth',1);
end
set([ax1; ax2],'Visible','off')
annotation('textbox',...
    [0.378333333333335,0.939130434782609,0.2,0.05],....
    'String','Neuron Traces','EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',15);
annotation('textbox',...
    [0.929166666666668,0.955460869565217,0.066666666666667,0.0391],....
    'String','Average Spike','EdgeColor','none','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',12);
tAx1 = axes('Position',[0.025 0.95-0.4/length(neurPlot)*A 0.90 0],'TickDir','out','XMinorTick','on');
linkaxes([ax1; tAx1])
xticklabels(xticks(tAx1)/acqRate)
xlabel(tAx1,'Time (s)')
%ax1(1).YLim = [-2.5e-4 2.5e-4]*2; % Use predetermined scale.
ax1(1).YLim = [-20 10]; % Use predetermined scale.

tAx2 = axes('Position',[0.94 0.95-0.4/length(neurPlot)*A 0.05 0],'TickDir','out','XMinorTick','on');
linkaxes([ax2; tAx2])
xticklabels(xticks(tAx2)/acqRate*1000)
xlabel(tAx2,'Time (ms)')
%ax2(1).YLim = [-2.5e-4 1e-4]; % Use predetermined scale.
ax2(1).YLim = [-20 10]; % Use predetermined scale.

set(hScOrig,'Visible','off')

%% Plot variance explained

ax3 = axes('Position',[0.05 0.05 0.3 0.3]);
plot(cumsum(expVar),'-o'); hold on
hVarThr = plot(get(gca,'xlim'),repelem(varThr,2),'--','LineWidth',2);
xlabel('Principal Component')
ylabel('Percent Variance')
ylim([0 100])
xlim([1 size(data,2)])
title('Cumulative Variance Explained')
legend(["Variance Explained" "Denoising Threshold"],'Location','best')

%% Plot covariance matrices

ax4(1) = axes('Position',[0.4 0.05 0.3 0.3]);
absCovOrig = padarray(abs(cov(data)),[1 1],'post');
covMatPlotOrig = surf(absCovOrig,'EdgeColor','none');
colorbar; axis image; view(0,90);
clim = caxis;
ylabel('Neuron #')
xlabel('Neuron #')
title('Covariance of Original Data')

ax4(2) = axes('Position',[0.70 0.05 0.3 0.3]);
absCovDenoised = padarray(abs(cov(denoised)),[1 1],'post');
covMatPlotDenoised = surf(absCovDenoised,'EdgeColor','none');
axis image; view(0,90);

%caxis(clim);
colorbar;
%ylabel('Neuron #')
yticks([])
xlabel('Neuron #')
title('Covariance of Denoised Data')

%% User interface panel
p2 = uipanel(f1,'Position',[0 0.4 1 0.1],'Title','Options');

% Checkboxes for toggling traces and spike times
cb1 = uicontrol(p2,'Style','checkbox','Position',[10 30 160 15],...
    'String','Show Original Traces','Value',1);
cb1.Callback = {@traceToggle,[hOrig hAvgOrig]};
cb2 = uicontrol(p2,'Style','checkbox','Position',[10 10 160 15],...
    'String','Show Denoised Traces','Value',1);
cb2.Callback = {@traceToggle,[hDenoised hAvgDenoised]};
cb3 = uicontrol(p2,'Style','checkbox','Position',[170 30 160 15],...
    'String','Show Original Spike Times','Value',0);
cb3.Callback = {@traceToggle,hScOrig};
cb4 = uicontrol(p2,'Style','checkbox','Position',[170 10 160 15],...
    'String','Show Denoised Spike Times','Value',1);
cb4.Callback = {@traceToggle,hScDenoised};

% Slider for neurons plotted
sld1 = uicontrol(p2,'Style','slider','Position',[110 50 220 15],...
    'Value',10,'Min',10,'Max',size(data,2),'SliderStep',[0.05 0.3]);
sld1.Callback = {@scrollTraces};
sldTx = uicontrol(p2,'Style','text','Position',[10 50 100 15],...
    'String','Neurons displayed');

% Editable variance explained threshold
ed1 = uicontrol(p2,'Style','edit','Position',[460 10 50 40],...
    'String',varThr,'FontSize',18);
ed1.Callback = {@chgVarThr,true};
ed1Tx1 = uicontrol(p2,'Style','text','Position',[430 50 100 15],...
    'String','Variance Threshold');

% Editable number of PCs removed threshold
ed2 = uicontrol(p2,'Style','edit','Position',[575 10 50 40],...
    'String',nPCs,'FontSize',18);
ed2.Callback = {@chgVarThr,false};
ed2Tx1 = uicontrol(p2,'Style','text','Position',[550 50 100 15],...
    'String','PCs Removed');

% Processing message
ed2Tx1 = uicontrol(p2,'Style','text','Position',[650 20 250 30],...
    'String','Denoising completed.','FontSize',15);
keyboard
%% Nested Functions
    function [denoised,nPCs] = removeNoise(nPCs)
        % This function removes noise by removing a given number of
        % Principal Components before reconstructing the data.
        if nargin~=0 % If a number of PCs is given as an input, remove that many PCs.
            eigZeroed = eigV;
            eigZeroed(:,1:nPCs) = 0;
            
            % Reconstruct centered data without removed components.
            denoised = Z * eigZeroed';
            
            % Complete reconstruction by adding back the means of each trace, that the
            % PCA algorithm subtracts when centering the data.
            denoised = denoised + repmat(mu,size(Z,1),1);
            denoised = zscore(denoised);
            
        else % If no number of PCs is given, find the optimal number to remove. This is done by minimizing the total covariance of the z-scored denoised data.
            
            zCovs = nan(size(data,2),1);
            zCovs(1) = sum(abs(cov(data)),'all'); % Total covariance of z-scored data.
            for A = 1:size(data,2)-1 % Iterate through removing all numbers of PCs, and determine resulting total covariance in each case.
                eigZeroed = eigV;
                eigZeroed(:,1:A) = 0;
                
                % Reconstruct centered data without removed components.
                denoised = Z * eigZeroed';
                
                % Complete reconstruction by adding back the means of each trace, that the
                % PCA algorithm subtracts when centering the data.
                denoised = denoised + repmat(mu,size(Z,1),1);
                denoised = zscore(denoised);
                
                zCovs(A+1) = sum(abs(cov(denoised)),'all'); % Total covariance of z-scored denoised data.
            end
            % Find minimum covariance and number of PCs to remove.
            [~,idx] = min(zCovs); 
            nPCs = idx - 1; % Need to subtract one because the first element of zCovs corresponds to removal of 0 PCs.
            
            % Recursive step. Call function again with number of PCs above
            % as an input.
            denoised = removeNoise(nPCs);
        end
    end

    function traceToggle(hObject,eventdata,handle)
        set(handle,'Visible',eventdata.Source.Value)
    end

    function scrollTraces(hObject,eventdata)
        sliderNeur = round(eventdata.Source.Value);
        neurPlot = sliderNeur-9:sliderNeur;
        updateTraces
    end

    function updateTraces
        for A = 1:length(neurPlot)
            set(neurLabel(A),'String',neurPlot(A))
            set(hOrig(A),'YData',data(:,neurPlot(A)))
            set(hDenoised(A),'YData',denoised(:,neurPlot(A)))
            
            tempSTO = orig_spiketimes(:,neurPlot(A));
            tempSTO(tempSTO==0) = [];
            tempSTO = tempSTO + spTimeOffset;
            set(hScOrig(A),'XData',tempSTO,'YData',data(tempSTO,neurPlot(A)))
            tempSTD = denoised_spiketimes(:,neurPlot(A));
            tempSTD(tempSTD==0) = [];
            tempSTD = tempSTD + spTimeOffset;
            set(hScDenoised(A),'XData',tempSTD,'YData',denoised(tempSTD,neurPlot(A)))
            
            set(hAvgOrig(A),'YData',orig_avgSp(:,neurPlot(A)))
            set(hAvgDenoised(A),'YData',denoised_avgSp(:,neurPlot(A)))
        end
    end

    function chgVarThr(hObject,eventdata,useVarNested)
        ed2Tx1.String = 'Processing...';
        drawnow
        if useVarNested
            varThr = str2double(eventdata.Source.String);
            nPCs = find(cumsum(expVar)>varThr,1);
        else
            nPCs = str2double(eventdata.Source.String);
            varThr = sum(expVar(1:nPCs));
        end
        
        [denoised,~] = removeNoise(nPCs);
        
        [denoised_spiketimes,denoised_spTraces,~,~]=detect_spikes_vsd(denoised);
        
        denoised_avgSp = cellfun(@(x) mean(x,1),denoised_spTraces,'uni',0);
        denoised_avgSp = cell2mat(denoised_avgSp')';
        
        
        updateVarThrPlot
        updateCovPlot
        updateTraces
        ed1.String = varThr;
        ed2.String = nPCs;
        ed2Tx1.String = 'Denoising completed.';
    end

    function testdisp(hObject,eventdata)
        disp(eventdata.Source.Value)
    end

    function updateVarThrPlot
        set(hVarThr,'YData',[varThr varThr])
    end

    function updateCovPlot
        absCovDenoised = padarray(abs(cov(denoised)),[1 1],'post');
        set(covMatPlotDenoised,'CData',absCovDenoised,'ZData',absCovDenoised)
    end
end