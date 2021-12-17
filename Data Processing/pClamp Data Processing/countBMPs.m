function [iBMPs,others] = countBMPs(thresholder)

% This function counts BMPs from the thresholder variable and outputs rates
% of iBMPs and other types of pattern per minute.

tWindow = 30; % Time window to consider for BMPs (in minutes).

isIBMP = thresholder.isIBMP(thresholder.protStart < tWindow*60);

iBMPs = sum(isIBMP)/tWindow;
others = sum(~isIBMP)/tWindow;