% This script is the main pipeline for NNMF analyses. 
% It requires the main dataset structure ("data").

%% Define variables for other sections.

% Folder dataset.
dPath = "E:\Renan\Operant Conditioning\blinded\16-Dec-2021Costa_et_al_complete_dataset.mat"; % Replace with local path.

% Load data.
load(dPath);

% Folder to save final dataset
saveFolder = 'E:\Renan\Operant Conditioning\blinded'; % Replace with local path.

%% Run NNMF

for ii = 1:numel(data)
    preSpT = data(ii).pretest.spikeTimes;
    postSpT = data(ii).posttest.spikeTimes;
    
    [NNMFpre,NNMFpost,~] = run_nnmf(preSpT,postSpT);
    
    data(ii).pretest.NNMF = NNMFpre;
    data(ii).posttest.NNMF = NNMFpost;
end

%%  Calculate correlation with BMP phases.

for ii = 1:numel(data)
    NNMFpre = data(ii).pretest.NNMF;
    NNMFpost = data(ii).posttest.NNMF;
    
    NNMFpre.protCorr = nnmfBMPCorr(data(ii).pretest.protStartEndTimes,NNMFpre.H);
    NNMFpre.retCorr = nnmfBMPCorr(data(ii).pretest.retStartEndTimes,NNMFpre.H);
    
    NNMFpost.protCorr = nnmfBMPCorr(data(ii).posttest.protStartEndTimes,NNMFpost.H);
    NNMFpost.retCorr = nnmfBMPCorr(data(ii).posttest.retStartEndTimes,NNMFpost.H);
    
    data(ii).pretest.NNMF = NNMFpre;
    data(ii).posttest.NNMF = NNMFpost;
end

%% Resample factors over a duration-normalized BMP

for ii = 1:numel(data)
    NNMFpre = data(ii).pretest.NNMF;
    NNMFpost = data(ii).posttest.NNMF;
    
    preProtFac = resampleFactors(data(ii).pretest.protStartEndTimes,NNMFpre.H);
    preRetFac = resampleFactors(data(ii).pretest.retStartEndTimes,NNMFpre.H);
    
    postProtFac = resampleFactors(data(ii).posttest.protStartEndTimes,NNMFpost.H);
    postRetFac = resampleFactors(data(ii).posttest.retStartEndTimes,NNMFpost.H);
    
    NNMFpre.H_BMP = cat(1,preProtFac,preRetFac);
    NNMFpost.H_BMP = cat(1,postProtFac,postRetFac);
    
    data(ii).pretest.NNMF = NNMFpre;
    data(ii).posttest.NNMF = NNMFpost; 
end

%% Match factors across preparations using pretest time courses.

for ii = 1:numel(data)
    avgH(:,:,:,ii) = mean(data(ii).pretest.NNMF.H_BMP,2); % Take mean of BMPs.
end

[facIds,~,facCorr] = matchFactors(avgH); % Obtain factor indexes

% Optional: visualize hierarchical clustering of factors.
if false
    factor_hClus(facCorr,size(avgH,3))
end
    

% Label factors
facLabels = ["retraction" "protraction"]; % if the number of factors were >2, more labels would be needed).
for ii = 1:numel(data)
    data(ii).pretest.NNMF.facLabels = facLabels(facIds(:,ii));
    data(ii).posttest.NNMF.facLabels = facLabels(facIds(:,ii));
end

%% Calculate phase(time) and magnitude of peak factor recruitment.

for ii = 1:numel(data)
    for testTag = ["pretest" "posttest"]
        avgH = mean(data(ii).(testTag).NNMF.H_BMP,2); % Take mean of BMPs.
        
        [mag,phase] = max(avgH,[],1);
        phase = phase./size(avgH,1); % Convert to normalized time scale.
        
        data(ii).(testTag).NNMF.peakMagnitude = squeeze(mag);
        data(ii).(testTag).NNMF.peakPhase = squeeze(phase);
    end
end

%% Optional: visualize contributions of each neuron to each factor.

if false
    facLabels = ["retraction" "protraction"];
    facColors = [85 136 198;... % Color for factor 1 (retraction)
    247 147 29]/256; % Color for factor 2 (protraction)
    for ii = 1:numel(data)
        for iii = 1:2
            facIdx = strcmp(data(ii).pretest.NNMF.facLabels,facLabels(iii));
            W = data(ii).pretest.NNMF.W(:,facIdx);
            order = data(ii).pretest.NNMF.cellOrderIdx;
            gImage = data(ii).posttest.gImage;
            det = data(ii).posttest.detKernels;
            color = facColors(iii,:);
            labelNeuronsNNMF(W,order,gImage,det,color);
        end
    end
end

%% Save data structure.

save(fullfile(saveFolder,strcat(date,'Costa_et_al_complete_dataset')),'data')