% This script is the main pipeline for data processing. It is provided here
% purely for illustration. The script cannot be run as it requires raw data
% and some functions may be missing.

%% Define variables for other sections. Always run this section first.

% Folder containing all, or any subset of, the experiments.
masterFolder = 'E:\Renan\Operant Conditioning\blinded\*';

% List .da files in experiment folders.
listf = dir(fullfile(masterFolder, '1*.da'));

% Load blinding key to assemble final dataset.
load('E:\Renan\Operant Conditioning\blinded\key.mat');

% Initialize final dataset data structure.
for A = key(:)'
    data(A).expLabel = keyExpLabels(key==A);
end

% Folder to save final dataset
saveFolder = 'E:\Renan\Operant Conditioning\blinded';

%% Find the optimal frame for drawing the kernels. This step will also create the .mat files on which data will be saved at other sections.

% Find optimal frame for each file.
for A = 1:length(listf)
    
    trial = listf(A).name(1:end-3); % Remove .da to obtain trial name.
    folder = listf(A).folder;
    try % In case there is an error on a specific file.
        find_kframe(folder,trial);
    catch
        display(['Failed to obtain optimal frame.' newline 'Folder: ' folder newline 'Trial: ' trial])
    end
end

display('Optimal frames for kernel drawing have been saved. Please, draw kernels for each file before proceeding.')

%% Extract raw data.

% Extract data for each file.
for A = 1:length(listf)
    
    trial = listf(A).name(1:end-3); % Remove .da to obtain trial name.
    folder = listf(A).folder;
    if isempty(dir(fullfile(folder,[trial '*.det']))) % Determine whether .det file is present.
        display(['No .det kernel file found for trial ' trial '. Data could not be extracted.'])
        continue
    else
        try % In case there is an error on a specific file.
            extract_data(folder,trial);
        catch
            display(['Failed to extract data.' newline 'Folder: ' folder newline 'Trial: ' trial])
        end
    end
end
%% Detect Spikes

% Detect spikes for each file.
for A = 1:length(listf)
    expID = floor((A+1)/2);
    if ~ismember(expID,key)
        continue
    end
    trial = listf(A).name(1:end-3); % Remove .da to obtain trial name.
    if strcmp(trial(3),'1')
        testTag = "pretest";
    elseif strcmp(trial(3),'2')
        testTag = "posttest";
    end
    folder = listf(A).folder;
    matF = dir(fullfile(folder,[trial 'pre*.mat']));
    try % In case there is an error on a specific file.
        
        load(fullfile(folder,matF.name),'vsddata','frame_pic','final_det')
        vsdFiltered = vsd_ellip(vsddata);
        vsdDenoise = pca_denoise(vsdFiltered);
        
        [spikeTimes,~,~,~]=detect_spikes_vsd(vsdDenoise);
        data(expID).(testTag).spikeTimes = spikeTimes;
        data(expID).(testTag).gImage = frame_pic;
        data(expID).(testTag).detKernels = final_det;
        
        clearvars('spikeTimes','frame_pic','final_det')
    catch
        display(['Failed to detect spikes.' newline 'Folder: ' folder newline 'Trial: ' trial])
    end
    
end

%% Obtain BMP times during imaging window.
% Note: imaging data is acquired using Neuroplex. Nerve recordings are
% acquired using both Neuroplex and pClamp, but pClamp data are used
% due to longer duration and higher acquisition rate, as well as for 
% consistency.
for A = 29:length(listf)
    expID = floor((A+1)/2);
    if ~ismember(expID,key)
        continue
    end
    trial = listf(A).name(1:end-3); % Remove .da to obtain trial name.
    if strcmp(trial(3),'1')
        testTag = "pretest";
    elseif strcmp(trial(3),'2')
        testTag = "posttest";
    end
    folder = listf(A).folder;
    thrF = dir(fullfile(folder,strcat(testTag,"thresholder.mat")));
    thr = load(fullfile(folder,thrF.name));
    
    [~, ~, ~, ~, ~, acqTimes] = findAcqTimes(folder,strcat(testTag,".abf"),strcat(trial,".da"),thr.recStart(end),thr.recEnd(end));
    
    [protStartEndTimes,retStartEndTimes,~,~] = findBMPTimes(acqTimes,thr);
    
    data(expID).(testTag).protStartEndTimes = protStartEndTimes;
    data(expID).(testTag).retStartEndTimes = retStartEndTimes;   
    
end


%% Count BMPs over full 30 min pretest and posttest nerve recordings (acquired using pClamp).

for A = 1:length(listf)
    expID = floor((A+1)/2);
    if ~ismember(expID,key)
        continue
    end
    trial = listf(A).name(1:end-3); % Remove .da to obtain trial name.
    if strcmp(trial(3),'1')
        testTag = "pretest";
    elseif strcmp(trial(3),'2')
        testTag = "posttest";
    end
    folder = listf(A).folder;
    thrF = dir(fullfile(folder,strcat(testTag,"thresholder.mat")));
    thr = load(fullfile(folder,thrF.name));
    
    [rateIBMPs,rateOtherBMPs] = countBMPs(thr);
    
    data(expID).(testTag).rateIBMPs = rateIBMPs;
    data(expID).(testTag).rateOtherBMPs = rateOtherBMPs;
    
    clearvars('thr','rateIBMPs','rateOtherBMPs')
end

%% Organize data structure.

emptyIdx = cellfun(@isempty,{data.expLabel});
data(emptyIdx) = [];
[~,sortIdx] = sort([data.expLabel]);
data = data(sortIdx);

%% Save data structure.

save(fullfile(saveFolder,strcat(date,'Costa_et_al_complete_dataset')),'data')