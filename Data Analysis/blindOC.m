%% Parameters
folder = 'E:\Renan\Operant Conditioning';

exp = {'18-04-17','18-04-20';'18-10-15','18-10-16';'18-10-22','18-10-24';'18-11-09','18-11-12';'18-11-26','18-12-03';'19-02-22','19-02-26';'19-06-13','19-06-18'};

dest = 'blinded';

%% Fail-safe checks
expisfolder = cellfun(@isfolder,fullfile(folder,exp));
if isfolder(fullfile(folder,dest))
    error('Aborted: Destination folder already exists.')
elseif ~all(expisfolder(:))
    error('Aborted: At least one experiment folder does not exist.')
elseif isfile(fullfile(folder,[datestr(date,'yy-mm-dd') '_Blinding_Pairing_Key.txt'])) ...
        || isfile(fullfile(folder,[datestr(date,'yy-mm-dd') '_Blinding_Full_Key.txt']))
    error('Aborted: At least one Blinding Key already exists for the current date in the specified folder.')
end

%% Move experiments from origin to blinded folders
blindFnames = num2str([randperm(numel(exp))]','%02d');
blindFnames = string(blindFnames);
blindFnames = reshape(blindFnames,size(exp));

for A = 1:numel(exp)
    movefile(fullfile(folder,exp{A}),fullfile(folder,dest,blindFnames{A}))
end
%% Generate complete key and save
fullKey = strcat(blindFnames(:),repmat(" is experiment ",numel(exp),1),exp(:),repmat(" \n ",numel(exp),1));
fullKey = strjoin(fullKey);

a=fopen(fullfile(folder,[datestr(date,'yy-mm-dd') '_Blinding_Full_Key' ...
    '.txt']) ...
    ,'w');
fprintf(a,fullKey);

%% Generate key that reveals which experiment is paired to which, but maintains contingent and yoked blinding.
sortedBlindFnames = sort(blindFnames,2); % Sorting columns then rows to ensure the order is not revealing.
sortedBlindFnames = sortrows(sortedBlindFnames); % This needs to be 'sortrows' instead of 'sort'. 'sort' breaks the pairings.
pairsKey = strcat(sortedBlindFnames(:,1),repmat(" is paired to ",size(sortedBlindFnames,1),1),sortedBlindFnames(:,2),repmat(" \n ",size(sortedBlindFnames,1),1));
pairsKey = strjoin(pairsKey);

b=fopen(fullfile(folder,[datestr(date,'yy-mm-dd') '_Blinding_Pairing_Key' ...
    '.txt']) ...
    ,'w');
fprintf(b,pairsKey);

fclose all;

disp('Experiments succesfully blinded.')
disp(['Blinded experiments location: ' convertStringsToChars(fullfile(folder,dest))])
disp(['Blinding keys location: ' folder])