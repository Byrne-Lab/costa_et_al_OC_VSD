function rsF = resampleFactors(fIdx,fMetric,rS)
% Resamples factors and outputs as time x BMPs x factors array.
% "fIdx"    One of the outputs of findBMPTimes. (I.e., indexes for one
%           phase). Can be either logical (Total time points x BMPs), or
%           linear (First row: phase start, second row: phase end, columns:
%           BMPs).
% fMetric   is the H output of seqNMF (i.e., time course of a factor).
% rS        Number of data points to use for each phase of each BMP.

if nargin<3
    rS = 180;
end
downSampFac = 10; % Must be the same as used for NNMF.

% Convert fIdx to logical if it is linear.
if ~islogical(fIdx)
    fIdxLin = fIdx;
    fIdx = false(size(fMetric,2)*downSampFac,size(fIdxLin,2));
    for ii = 1:size(fIdxLin,2)
        fIdx(fIdxLin(1,ii):fIdxLin(2,ii),ii) = true;
    end
end

fIdx(:,~any(fIdx,1)) = [];
fIdx = downsample(fIdx,downSampFac);
for ii = 1:size(fIdx,2)
    for kk = 1:size(fMetric,1)
        rawMetric = fMetric(kk,fIdx(:,ii));
        rsF(:,ii,kk) = interpRF(rawMetric,rS,length(rawMetric)); % Changed from resample to interpRF. '25-Jan-2021'
    end
end
end


function y = interpRF(x,n,d) % Simple interpolation from rational fraction of sampling interval. Same inputs as the MATLAB resample function, but performs only simple linear interpolation (i.e, no low pass filter).
y = interp1(...
    linspace(0,1,length(x))...
    ,x...
    ,linspace(0,1,length(x)*n/d)...
    );
end