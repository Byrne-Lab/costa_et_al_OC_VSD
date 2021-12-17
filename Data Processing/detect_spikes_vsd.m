function [spTimes,spTrace,spLogical,detect] = detect_spikes_vsd(data)

%% Parameters

acqRate = 1;        % Acquisition rate used for the recording, in kHz.
W = [15,60];        % Window before and after peak for visualize spike.
ra = 10;            % Duration (ms) to re-arm the spike detection.
thr = 5;            % Amplitude threshold for spike detection. In standard deviations of the estimated noise (details in code).

%%

tRange = -W(1):W(2);

% Detect Spikes
spLogical = detection(data);

% Build spike times output
spTimes = getSpTime(spLogical);

% Build spike traces output (cell array of individual spike traces for each
% neuron)
spTrace = getSpTrace(spLogical);

% Temporarily included for back compatibility. Useless output.
detect = [];

%% Nested functions
    function detect_all = detection(data)
        dataz = zscore(data);
        detect_all = false(size(dataz));
        for A = 1:size(dataz,2)
            %             if ~noiseType
            %             noiseEst = std(dataz(:,A));
            %             else
            noiseEst = median(abs(dataz(:,A))/0.6745); % Estimation of
            % the noise level. Equation obtained from Quiroga RQ, Nadasdy
            % Z, Ben-Shaul Y (2004) Unsupervised Spike Detection and
            % Sorting with Wavelets and Superparamagnetic Clustering.
            % Neural Comput 16:1661–1687 Available at:
            % http://www.mitpressjournals.org/doi/10.1162/089976604774201631.
            % Equation originally introduced by Donoho DL, Johnstone IM (1994) Ideal
            % spatial adaptation by wavelet shrinkage. Biometrika
            % 81:425–455 Available at:
            % https://academic.oup.com/biomet/article/81/3/425/256924.
            %             end
            warning('off','signal:findpeaks:largeMinPeakHeight');
            
            [~,spTemp] = findpeaks(-dataz(:,A),'MinPeakHeight',thr*noiseEst,'MinPeakDistance',round(ra*acqRate));
            % Note: For some reason, if the firing rate is given to
            % findpeaks, its performance is worse. Might be an offset
            % issue due to how findpeaks generates the time vector "x".
            % Thus, "ra" is converted to datapoints directly and then
            % used as an input.
            
            detect_all(spTemp,A) = true;
        end
    end

    function spTimes = getSpTime(spBinary)
        % Build spike times output
        spTimes = zeros(max(sum(spBinary,1)),size(data,2));
        for A = 1:size(data,2)
            spTimes(1:sum(spBinary(:,A)),A) = find(spBinary(:,A));
        end
    end

    function spTrace = getSpTrace(spBinary)
        % Build spike traces output (cell array of individual spike traces for each
        % neuron)
        spTrace = cell(1,size(data,2));
        for A = 1:size(data,2)
            spTrace{A} = nan(sum(spBinary(:,A)),length(tRange));
            for B = 1:sum(spBinary(:,A))
                tWin = spTimes(B,A)+tRange;
                tWin(tWin<1) = 1; tWin(tWin>size(data,1)) = size(data,1); % Trim window beyond data limits. For plotting consistency, points beyond limits are made equal to the first or last point.
                spTrace{A}(B,:) = data(tWin,A);
            end
        end
    end

end