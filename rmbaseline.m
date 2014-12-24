
function [out,bl] = rmbaseline(bx,usepts,varargin)


% Fits a polynomial to the spectrum and normalizes by the fitted value.

% C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

smoothing_method = 'polynomial';

if nargin < 2 || isempty(usepts)
    usepts = true;
end
if nargin < 3
    polyord = 8;
elseif isnumeric(varargin{1})
    polyord = varargin{1};
else
    
    smoothing_method = varargin{1};
    if nargin < 4,
        smoothbw = 10;
    else
        smoothbw = varargin{2};
    end
end

mn = mean(abs(bx.blrep),1);
kp = mn~=0 & usepts;
switch smoothing_method
    case 'polynomial'
        p = polyfit(find(kp)./length(mn),log(mn(kp)),polyord);
        bl = polyval(p,(1:length(mn))/length(mn));
    case 'moving_average'
        g = rectwin(smoothbw)';
        bl = log(convn(mn.*kp,g,'same')./convn(kp,g,'same'));
        
end



out = bx.blrep./repmat(exp(bl),size(bx.blrep,1),1);

