%% Function for finding times of frame acquisition. This function returns the clampfit indexes for each of the frames recorded in neuroplex.
% stimIdx are the indexes of the first data point acquired by Neuroplex
% after each stimulus.
% NplxTimes are the times at which each data point in each BNC channel was acquired by Neuroplex
% relative to Clampfit time. stimTimes are the times at which each stimulus
% occurs in Clampfit time. stimDur is the median stimulus duration and
% stimInt is the median stimulus intensity. NplxAvgTime are estimated for
% acquisition of each frame by Neuroplex, based on the median of all BNC
% channels.
function [stimIdx, NplxTimes, stimTimes, stimDur, stimInt, NplxAvgTime] = findAcqTimes(folder,CfitFn,NplxFn,recStart,recEnd)
%% Parameters

% BNC channel carrying the sync wave on Neuroplex.
NplxSyncBNC = 1;

% BNC channels carrying nerve traces on each file
NplxNerves = 2:8;
CfitNerves = 3:9;

% Downsampling factor for Clampfit to speed up the process. (Probably not a
% good idea, to be removed. Set to 1 for no downsampling.)
downSampFac = 1;

%% Open Neuroplex file

fid = fopen(fullfile(folder, NplxFn), 'r', 'l');

%% Read header info

% Number of columns
fseek (fid, 2*384, 'bof');
par.dim(1) = double(fread(fid,1,'*int16'));

% Number of rows
fseek (fid, 2*385, 'bof');
par.dim(2) = double(fread(fid,1,'*int16'));

% Number of frames (time points)
fseek (fid, 2*2000, 'bof');
par.dim(3) = sum(double(fread(fid,2,'*int16')) .* [1; 32000]); % Don't know why, but this works for large data files. (Maybe write directly to disk large files? Not sure.)

% Frame interval (ms)
fseek (fid, 2*388, 'bof');
par.framedt = double(fread(fid,1,'*int16')) / 1000;
fseek (fid, 2*390, 'bof');
dividingFactor = double(fread(fid,1,'*int16'));
if par.framedt >= 10
    par.framedt = par.framedt * dividingFactor;
end
par.framedt = par.framedt / 1000; % Convert from ms to s

% Aquisition ratio (electrical over optical)
fseek (fid, 2*391, 'bof');
par.acqRatio = double(fread(fid,1,'*int16'));
if par.acqRatio == 0
    par.acqRatio = 1;
end

%% Load Neuroplex data

fseek(fid, 5120+2*prod(par.dim), 'bof');      % Move to the beginning of BNC data in the data file
NplxBNC = fread(fid, par.dim(3)*par.acqRatio*8, 'int16');  % read all BNC data regardless of num_frame (to be simple) as a single vector
NplxBNC = reshape(NplxBNC, par.dim(3)*par.acqRatio, 8);    % split to 8 channels

% Normalize
NplxBNC = normBNC(NplxBNC);

% Load first 1000 frames and compute average light intensity. (To obtain
% time of shutter opening.)
fseek(fid, 5120, 'bof'); 
NplxLI = fread(fid, 1000*par.acqRatio*prod(par.dim(1:2)), 'int16');
NplxLI = mean(reshape(NplxLI,par.dim(1)*par.dim(2),1000),1);

fclose('all');
%% Load Clampfit data
[CfitBNC,CfitFramedt,h]=abfload(fullfile(folder, CfitFn),'doDispInfo',0);

% Trim Clampfit data to recorded window
if nargin < 5
    recStart = 0;
    recEnd = size(CfitBNC,1)*CfitFramedt/1e6-CfitFramedt/1e6;
else
    recStart = recStart - 1; % -1 second offset to capture shutter opening.
end
recStartIdx = round(recStart/CfitFramedt*1e6 + 1); 
recEndIdx = round(recEnd/CfitFramedt*1e6 + 1);
CfitBNC = CfitBNC(recStartIdx:recEndIdx,:);

% Downsample trace for speed.
CfitBNC = downsample(CfitBNC,downSampFac);
CfitFramedt = CfitFramedt*downSampFac;

% Reshape data to concatenate trials.
CfitBNC = permute(CfitBNC, [1 3 2]);
CfitBNC = reshape(CfitBNC,[],size(CfitBNC,3),1);

% Find box 1 stimulus trace index.
box1Idx = cellfun(@(x) strcmp('STIMBOX1',x),h.recChNames);

CfitSyncBNC = cellfun(@(x) strcmp('SYNC',x),h.recChNames);
CfitSyncBNC = find(CfitSyncBNC);
CfitClockBNC = cellfun(@(x) strcmp('Clock Trig',x),h.recChNames);
CfitShutBNC = cellfun(@(x) strcmp('Shutter',x),h.recChNames);

% Full resolution stimulus trace.
stimTrace = CfitBNC(:,box1Idx);

% Normalize.
CfitBNC = normBNC(CfitBNC);

% Time vector (s)
CfitFramedt = CfitFramedt/1e6; % Convert from us to s.
CfitTime = 0:CfitFramedt:size(CfitBNC,1)*CfitFramedt-CfitFramedt;

%% Find acquisition times

% Index clock pulses from the camera.
clockBNC = CfitBNC(:,CfitClockBNC);
[~,clockPeaks] = findpeaks(diff(clockBNC),'MinPeakHeight',0.1);
%CfitIdx = CfitIdx + 1; % Correcting for "diff"

% Obtain discrepancy constant between the two clocks (Clampex over
% Neuroplex; i.e., how much clampex time is equivalent to "1" neuroplex
% time -- don't let the equation confuse you)
clkDisc = mean(diff(clockPeaks))*CfitFramedt/par.framedt;

% Obtain rational fraction of the discrepancy (in data points, for
% resampling.
[nRF,dRF] = rat(clkDisc*par.framedt/CfitFramedt);

% Resample Neuroplex trace to Clampex data points.
%resNplxBNC = resample(NplxBNC,nRF,dRF); % Using interpolation instead of
%resampling. '02-Jan-2021'
resNplxBNC = interpRF(NplxBNC,nRF,dRF);

% Resample average light intensity.
%resNplxLI = resample(NplxLI,nRF,dRF); % Using interpolation instead of
%resampling. '02-Jan-2021'
resNplxLI = interpRF(NplxLI,nRF,dRF);

% Find times of shutter opening in Neuroplex (from light intensity) and in
% clampfit (shutter BNC output).
shutNplx = find(resNplxLI>1000,1); % Time when the average absolute light intensity is above 1000.
shutCfit = find(CfitBNC(:,CfitShutBNC)>-5,1);

% Find the estimated overall delay between recordings (Clampex to Neuroplex)
clkDelayCfitEst = shutCfit - shutNplx;
clkDelayCfitEst = clkDelayCfitEst + 45; % Add offset. 45 datapoints was determined empirically from the data. This is to account for the time it takes for the mechanical shutter to start opening.


% Find the delay between recordings for each channel (Clampex to Neuroplex)
[clkDelayCfit, gapSizeCfit, gapPositionCfit] = nestedFindDelay(resNplxBNC(:,[NplxSyncBNC NplxNerves])...
     ,CfitBNC(:,[CfitSyncBNC CfitNerves]),...
     clkDelayCfitEst);
clkDelay = clkDelayCfit*CfitFramedt/par.framedt; % clkDelayCfit is in Clampex data points, clkDelay in neuroplex datapoints.
gapSize = gapSizeCfit*CfitFramedt/par.framedt;
gapPosition = gapPositionCfit*CfitFramedt/par.framedt;
%clkDelay = median(clkDelay); % Set final delay as median of all channels.
%clkDelay = clkDelay(1); % Or just use the first. Mainly for testing, comment either this or previous line out.


% Calculate Neuroplex time base adjusted to Clampex data points. (For each
% nerve. Delays are different for each channel due to the time it takes to
% acquire and store data, and the order in which this is done by Neuroplex
% vs Clampfit.)

NplxTimes = repmat(0:size(NplxBNC,1)-1,length(clkDelay),1);
NplxTimes = (NplxTimes * clkDisc) + repmat(clkDelay',1,size(NplxTimes,2));
for A = 1:size(gapSize,2)
    if ~isnan(gapSize(A)) && ~isnan(gapPosition(A))
        gapIdx = NplxTimes(A,:)>gapPosition(A);
        NplxTimes(A,gapIdx) = NplxTimes(A,gapIdx) + gapSize(A);
    end
end
NplxTimes = NplxTimes * par.framedt; % Convert from Nplx data points (ms at 1 kHz) to s.

NplxTimes = NplxTimes + recStart; % Add offset from removed data.

NplxAvgTime = median(NplxTimes,1,'omitnan'); % Approximate overall time of acquisition of each frame.

%% Calculate times, and mean duration and intensity of stimuli

% Find onset of each stimulus
[~,stimStartIdx] = findpeaks(diff(stimTrace),'MinPeakHeight',1);
stimStartIdx = stimStartIdx + 1; % To correct for diff

if isempty(stimStartIdx) || ~any(stimTrace>1)  % To skip subsequent code if no stimulus is present.
    stimIdx = [];
    stimTimes = [];
    stimDur = [];
    stimInt = [];
else
    
    % Find offset of each stimulus
    [~,stimEndIdx] = findpeaks(diff(-stimTrace),'MinPeakHeight',1);
    stimEndIdx = stimEndIdx + 1; % To correct for diff
    
    % Stimulation times are the onsets of each stimulus on Clampfit Time
    stimTimes = CfitTime(stimStartIdx);
    
    % Calculate mean stimulus duration
    stimDur = CfitTime(stimEndIdx) - stimTimes;
    stimDur = mean(stimDur);
    
    % Calculate mean stimulus intensity
    stimAllIdx = [];
    for A = 1:length(stimStartIdx)
        stimAllIdx = cat(2,stimAllIdx,...
            stimStartIdx(A):stimEndIdx(A)-1);
    end
    stimInt = stimTrace(stimAllIdx);
    stimInt = mean(stimInt);
    
    % This is a seemingly convoluted but computationally inexpensive solution
    % to a simple problem. "stimTimes" is an array of the times of stimuli in the
    % clampfit trace (i.e., with the clampfit acq rate), but needs to be
    % converted to the equivalent index in the neuroplex traces. "NplxAvgTime"
    % is the time of each data point acquired by Neuroplex in the clampfit
    % trace (i.e., with the clampfit acq rate). Thus, the correct way to do
    % this is to find the values of
    % "NplxAvgTime" that are either identical to or the immediately next larger
    % value than each "stimTimes" value. The procedure below avoids the use of a
    % for loop for performing such comparisons.
    stimIdx = find(ismember(sort([stimTimes, NplxAvgTime]),stimTimes)); % "stimTimes"  and "NplxAvgTime"
    % are concatenated and sorted. The indexes in this new vector of the
    % values that are contained in the original "stimTimes" are the converted
    % "stimIdx" values to index the neuroplex trace, but a few corrections
    % are necessary.
    stimIdx = stimIdx([true, diff(stimIdx)>1]); % If the index difference is equal
    % to one, that means that the "stimIdx" value was repeated (i.e., already
    % contained by "NplxAvgTime" and thus needs to be removed.
    stimIdx = stimIdx - (0:(length(stimIdx)-1)); % Concatenating and sorting the
    % vectors pushes all values back by the number of values that have been
    % added until that point. Thus, it is necessary to subtract a vector of
    % this displacement to correct the final indexes.
end

%% Normalizing function. % Usage currently commented out. May or may not help. 03/03/2020
    function nBNC = normBNC(BNC)
        
        % Subtract mean.
        nBNC = BNC - repmat(mean(BNC),size(BNC,1),1);
        
        % Divide by range.
        nBNC = nBNC ./ (repmat(std(BNC),size(BNC,1),1));
    end

%%
    function [delay, gapSz, gapPos] = nestedFindDelay(x,y,delayEstimate) % Custom find delay function
        numPeaks = 2; % Fixed at 2 for now.
        
        %rho = nan(max(size(x,1),size(y,1),size(x,2)));
        [peaks, delay] = deal(nan(numPeaks,size(x,2)));
        for ii = 1:size(x,2)
            if all(isnan(x(:,ii))) || all(isnan(y(:,ii)))
                continue
            end
            rho = xcorr(x(:,ii),y(:,ii));
            [pks, locs] = findpeaks(rho,'SortStr','descend','NPeaks',numPeaks+1); 
            
            
            if ~isempty(pks)&&~isempty(locs)
                % Keep correlations that stand out only.
                pksKeep = -diff(pks);
                pksKeep = pksKeep > pksKeep(1)*0.1;
                
                pks = pks(1:end-1); % To remove extra peak.
                locs = locs(1:end-1); 
                
                pks(~pksKeep) = nan;
                locs(~pksKeep) = nan;
                
                %locs(pks<0.5*pks(1)) = nan; Commented out, using different
                %solution (above). '06-Jan-2021'
                peaks(:,ii) = pks;
                delay(:,ii) = locs;
            end
        end
        
        delay = -(delay - max(size(x,1),size(y,1))); % Offset. The zero lag position is at the center of the cross-correlation, so the largest vector length must be subtracted to obtain the lag. Signal flipped for consistency with output of the built-in finddelay function.
        
        % Reorder delays by placing delays that are close to the shutter
        % estimate on top.
        [~,reorderIdx] = sort(abs(delay-delayEstimate),1,'ascend');
        reorderIdx = sub2ind(size(reorderIdx),reorderIdx,repmat(1:size(reorderIdx,2),size(reorderIdx,1),1)); % Convert from row index to linear index.
        delay = delay(reorderIdx);
        
        % Use delay estimate if no peaks are close to shutter opening.
        if ~any(abs(delay-delayEstimate)<45)
            delay(~isnan(delay)) = delayEstimate;
            delay(2,:) = nan;
            warning('Using shutter opening to estimate delay between Clampfit and Neuroplex.')
        end
        
        gapPos = nan(1,size(x,2));
        if sum(~isnan(delay(2,:)))>5 % If at least 6 nerves show a second peak correlation, run sliding gap function.
            gapPos = nan(1,size(x,2));
            for ii = 1:size(x,2)
                tic;
                rhoGap = slidingGap(x(:,ii),y(:,ii),delay(:,ii));
                toc
                [~,gapPos(ii)] = max(rhoGap);
                
            end
        end
        gapSz = abs(diff(delay));
        gapPos(isnan(gapSz)) = nan;
        
        delay(2:end,:) = [];
        
    end

%% Function to compute correlation between vectors x and y following addition of a sliding gap to y. 
    function rho = slidingGap(x,y,delays)
        % delays are the possible delays obtained from the custom find delay function
        % stepSz is the number of steps to slide per iteration.
        % rho is set to have the same number of elements as x anx y
        % regardless of gap size.
        
        
        if all(isnan(delays(2:end)))
            rho = zeros(size(x));
            return
        end
        
        % gapSz is the size of the sliding gap 
        gapSz = abs(diff(delays));
        if length(gapSz)>1
            error('Multiple sliding gaps not currently implemented.')
        end
        
        % Align vectors based on first delay
        if delays(1)>0
            x = cat(1,zeros(delays(1),1),x);
        elseif delays(1)<0
            y = cat(1,zeros(-delays(1),1),y);
        end
        
        % Match vector lengths by appending zeros at end of smaller vector
        lengDiff = size(x,1)-size(y,1);
        if lengDiff>0
            y = cat(1,y,zeros(lengDiff,1));
        elseif lengDiff<0
            x = cat(1,x,zeros(-lengDiff,1));
        end
        
        if size(x,2)>1 || size(y,2)>1
            error('Sliding gap correlation can only be computed for vectors.')
        end
        
        gap = zeros(gapSz,1);
        yii = cat(1,y,gap);
        rho = nan(size(x));
        
        % Initiate variables for optimization
        stepSz = round(size(x,1)/50);
        pos = [1; size(x,1)-gapSz];
        chunkSz = abs(diff(pos));
        iter = 0;

        while stepSz>0 && iter<100
            for ii = pos(1):stepSz:pos(2)
                xii = cat(1,x(1:ii-1),gap,x(ii:end));
                rho(ii:ii+stepSz-1) = corr(xii,yii);
            end
            idx = find(rho==max(rho));
                   
            chunkSz = chunkSz/4;
            pos(1) = round(max([1 idx(1)-chunkSz/2]));
            pos(2) = round(min([size(x,1)-gapSz idx(end)+chunkSz/2]));
            stepSz = round(stepSz/4);
            iter = iter + 1;
            %disp(['Iteration: ' num2str(iter)])
        end
        
    end

    function y = interpRF(x,n,d) % Simple interpolation from rational fraction of sampling interval. Same inputs as the MATLAB resample function, but performs only simple linear interpolation (i.e, no low pass filter).
        y = interp1(...
                linspace(0,1,length(x))...
                ,x...
                ,linspace(0,1,length(x)*n/d)...
            );
    end
end