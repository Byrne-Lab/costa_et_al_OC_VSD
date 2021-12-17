function [protIdxLin,retIdxLin,protIdxLogic,retIdxLogic] = findBMPTimes(acqTimes,thresholder)
% This function outputs indexes for the protraction and retraction phases
% relative to the imaging window. It compares the BMP times during pClamp
% nerve recordings with the acquisition times for each datapoint recorded
% by Neuroplex during the imaging window.
%
% Variables:
% acqTimes      Last output of findAcqTimes ("NplxAvgTime"). Time of the
%               acquisition of each Neuroplex frame relative to the pClamp
%               recording.
% thresholder   Output of the BMP_thresholder script. Structure that 
%               contains start and end times for each phase of each BMP
%               captured by pClamp.


% Exclude BMPs that begin or end outside of the imaging window. Added
% on  '04-Feb-2021'
removIdx = (thresholder.protStart<acqTimes(1) | thresholder.retEnd>acqTimes(end));
thresholder.protStart(removIdx) = [];
thresholder.protEnd(removIdx) = [];
thresholder.retEnd(removIdx) = [];

[protIdxLogic,retIdxLogic] = ...
    deal(false(size(acqTimes,2),length(thresholder.protStart)));
for jj = 1:length(thresholder.protStart)
    protIdxLogic(:,jj) = acqTimes>=thresholder.protStart(jj) & acqTimes<=thresholder.protEnd(jj)';
    retIdxLogic(:,jj) = acqTimes>thresholder.protEnd(jj) & acqTimes<=thresholder.retEnd(jj)';
end

[protIdxLin,retIdxLin] = ...
    deal(nan(2,length(thresholder.protStart)));
for jj = 1:size(protIdxLogic,2)
    protIdxLin(1,jj) = find(protIdxLogic(:,jj),1,'first');
    protIdxLin(2,jj) = find(protIdxLogic(:,jj),1,'last');
    retIdxLin(1,jj) = find(retIdxLogic(:,jj),1,'first');
    retIdxLin(2,jj) = find(retIdxLogic(:,jj),1,'last');
end

