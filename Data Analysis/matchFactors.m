function [facIds,matches,rho] = matchFactors(factors)
% Function for matching factors across experiments on a normalized
% timescale.

% factors   is a four-dimensional array of normalized timepoints x
%           BMPs x factors x preparations. This variable should be an average
%           of BMPs, so the size of the second dimension should always
%           equal 1.
%

nFacs = size(factors,3);

factors = reshape(factors,size(factors,1),numel(factors)/size(factors,1)); 
% Reduce to two-dimensional array. Rows are timepoints, columns are factors*experiments. 
% (i.e., fac1exp1 fac2exp1 fac1exp2 fac2exp2, etc.)

rho = corr(factors); % Correlation is equivalent to centered cosine similarity.
%rho = factors'*factors./(vecnorm(factors)'*vecnorm(factors));% Cosine similarity (uncentered).
rhoIter = rho;

% Creating a mask to prevent within-preparation matches.
mask = false(size(rhoIter));
mask(1:size(mask,1)*nFacs+nFacs:end) = true;
mask = logical(conv2(mask,true(nFacs)));
mask = mask(1:size(rhoIter,1),1:size(rhoIter,2));

matches = false(size(rhoIter));

rhoIter(mask) = nan;
ii = 0;
while any(~isnan(rhoIter),'all')
    ii = ii+1;
    [~,idx] = max(rhoIter,[],'all','linear');
    matches(idx) = true;
    
    [row, col] = ind2sub(size(rhoIter),idx);
    
    posC = mod(col+nFacs-1,nFacs)+1; % To get the position of the matched factor relative to other factors.
    idxFacsC = 1:nFacs;
    idxFacsC = idxFacsC - posC;
    
    
    posR = mod(row+nFacs-1,nFacs)+1; % To get the position of the matched factor relative to other factors.
    idxFacsR = 1:nFacs;
    idxFacsR = idxFacsR - posR;
    
    rhoIter(row,col+idxFacsC) = nan;
    rhoIter(row+idxFacsR,col) = nan;
    rhoIter(col+idxFacsC,row) = nan;
    rhoIter(col,row+idxFacsR) = nan;% Also removing the symmetric (transposed) position.
    
    if ii>=size(rhoIter,1)*size(rhoIter,2)/nFacs
        warning('Iterations on while loop went over maximum expected. Potentially stuck.')
    end
end

% Generate IDs for indexing factors across preparations.
facIds = matches | matches' | eye(size(matches));
facIds = double(facIds(1:nFacs,:));

for ii = 1:nFacs
    facIds(ii,:) = facIds(ii,:).*ii; 
end

facIds = sum(facIds,1);
facIds = reshape(facIds,nFacs,size(facIds,2)/nFacs);

end

