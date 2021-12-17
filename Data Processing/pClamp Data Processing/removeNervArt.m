%% This function removes stimulation artifacts from nerve recordigns by setting the values of a window around the stimuli to 0.
function [cleanData, si, h] = removeNervArt(abfFile)

%% Parameters
%abfFile = 'C:\Data\Renan\Operant Conditioning\18-01-26\101training102.abf';

% Windows to erase in ms
preStimWind = 5; % Window preceding each stimulus.
box1Wind = 25; % Window following stimuli from stimulus Box 1.
box2Wind = 50; % Window following stimuli from stimulus Box 2.

%%

[data, si, h] = abfload(abfFile,'doDispInfo',0);

box1Idx = find(cellfun(@(x) strcmp('STIMBOX1',x),h.recChNames));
box2Idx = find(cellfun(@(x) strcmp('STIMBOX2',x),h.recChNames));


if any(data(:,box1Idx) > 1) || any(data(:,box2Idx) > 1)

% OLD APPROACH, now obsolete (2020-11-09)
%     box1Erase = data(:,box1Idx) > 1;
%     for A = 1:10:box1Wind/(si*10^-3)
%         box1Erase = box1Erase | circshift(box1Erase,10);
%     end
%     
%     box2Erase = data(:,box2Idx) > 1;
%     for A = 1:10:box2Wind/(si*10^-3)
%         box2Erase = box2Erase | circshift(box2Erase,10);
%     end
%     
%     for A = 1:10:preStimWind/(si*10^-3)
%         box1Erase = box1Erase | circshift(box1Erase,-10);
%         box2Erase = box2Erase | circshift(box2Erase,-1);
%     end
    
    box1Erase = convAndShift(data(:,box1Idx),box1Wind/(si*10^-3),preStimWind/(si*10^-3));
    box2Erase = convAndShift(data(:,box2Idx),box2Wind/(si*10^-3),preStimWind/(si*10^-3));

    box1Erase = repmat(box1Erase,1,size(data,2),size(data,3));
    box2Erase = repmat(box2Erase,1,size(data,2),size(data,3));
    
    box1Erase(:,[box2Idx,box1Idx],:) = false; % To prevent erasing the recorded stimulus itself.
    box2Erase(:,[box2Idx,box1Idx],:) = false; % To prevent erasing the recorded stimulus itself.
    
    cleanData = data;
    cleanData(box1Erase) = 0;
    cleanData(box2Erase) = 0;
else
    display('No stimulus in data')
    cleanData = data;
end

    function eraseIdx = convAndShift(boxVec,postWinLength,preWinLength)
        postWinLength = round(postWinLength);
        preWinLength = round(preWinLength);
        
        stimIdx = boxVec>1;
        winVec = true(postWinLength,1);
        
        eraseIdx = conv(stimIdx,winVec,'same');
        eraseIdx = eraseIdx~=0;
        
        pad = false(floor(postWinLength/2 - preWinLength),1);
        
        eraseIdx = cat(1,pad,eraseIdx(1:length(eraseIdx)-length(pad)));
    end

end
