function [phaseCorr] = nnmfBMPCorr(fIdx,fMetric)
% Calculates the correlations among each phase and each factor.
% "fIdx"    One of the outputs of findBMPTimes. (I.e., indexes for one
%           phase). Can be either logical (Total time points x BMPs), or
%           linear (First row: phase start, second row: phase end, columns:
%           BMPs).
% fMetric   is the H output of seqNMF (i.e., time course of a factor).

downSampFac = 10; % Must be the same as used for NNMF.


if ~islogical(fIdx)
    fIdxLin = fIdx;
    fIdx = false(size(fMetric,2)*downSampFac,size(fIdxLin,2));
    for ii = 1:size(fIdxLin,2)
        fIdx(fIdxLin(1,ii):fIdxLin(2,ii),ii) = true;
    end
end

phaseCorr = corr(fMetric',downsample(any(fIdx,2),downSampFac));
