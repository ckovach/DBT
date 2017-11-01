function x = iterz(x,thresh,side)

% xout = iterz(xin,[threshold=10],[side=0])
%
% Iterative z-score thresholding. At each iteration a z-score is computed
% on the columns of x, ignoring NaNs, and all values in x above the given z
% threshold are replaced with NaNs until no values in x are above the threshold. 
% If side = 0, thresholding is based on the magnitude of the z-score. Thresholding 
% is applied to positive values if side = 1 and negative values if side =
% -1.
%
%

% C Kovach 2017

if nargin < 2
    thresh = 10;
end
if nargin < 3
    side = 0;
end
    

if side ==0
    zscore = @(x)abs(x-repmat(nanmean(x),size(x,1),1))./repmat(nanstd(x),size(x,1),1);
elseif abs(side)==1;
    zscore = @(x)side*(x-repmat(nanmean(x),size(x,1),1))./repmat(nanstd(x),size(x,1),1);
else
    error('Side parameter must be 0 (two-sided), 1( high threshold ), or -1 (low threshold)')
end

z = zscore(x);

while(any(z(:)>thresh))
    x(z>thresh)=nan;
    z = zscore(x);
end
    
    