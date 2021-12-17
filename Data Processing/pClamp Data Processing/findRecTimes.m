function [recStart, recEnd, shutStart, shutEnd] = findRecTimes(abfFile)

%%Parameters

thr = 1; % 1 V change threshold.
%%
[data, si, h] = abfload(abfFile,'doDispInfo',0);

shutterIdx = find(cellfun(@(x) strcmp('Shutter',x),h.recChNames));

shutter = diff(data(:,shutterIdx));
shutter(circshift((shutter>thr)|(shutter<-thr),1)) = 0; % Remove adjacent true values, helps with instances where voltage change is not instantaneous relative to acquisition rate.

shutStart = find(shutter>thr)*si/10^6; % In seconds.
shutEnd = find(shutter<-thr)*si/10^6; % In seconds.

recStart = shutStart; % In order to get both shutter and recording times out.
recEnd = shutEnd;

% Remove short shutter openings (recordings should be ~2 min)
removIdx = ((recEnd - recStart) < 110) | ((recEnd - recStart) > 130);
recStart(removIdx) = [];
recEnd(removIdx) = [];

end

