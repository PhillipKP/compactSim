function [filterOrder] = hybrid18FiltFunc(S, numSpecFilt, specFiltersVis)


% The scatter matrix and the eigenvectors
G = S * S';
[V,~] = eig(G);

% flips the all principal components left to right, since the
V = fliplr(V);

for mInd = 1:numSpecFilt
    
    if mInd < 19
        
        % Select the first 6 spectral filters the are most parallel to the
        % first 6 endmembers
        
        % This is the principal component. But in reality we cannot implement
        % arbitrary spectral filters.
        fIdeal = V(:, mInd);
        
        % So we create a list of angles that each spectral filter makes with
        % the principal component vector for this measurment number
        cosAngList = computeCosine(fIdeal, specFiltersVis );
        
        % We then take the absolute value of the cos of the angles
        % Does the absolute value this matter? I'm not sure!
        absCosAngList = abs(cosAngList);
        
        
        % maxInd is the spectral filter that has the largest cosine with
        % the ideal spectral filter. Large cosine means their angles are likely
        % to be parallel
        [~, maxInd] = max(absCosAngList);
        
        % Add this index to the list.
        filterOrder(mInd) = maxInd;
        
    else
        
        hrInd = randNumNotInList(numSpecFilt, filterOrder);
        
        % Add this index to the list.
        filterOrder(mInd) = hrInd;
        
    end
    
    
end



end