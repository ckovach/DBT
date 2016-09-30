
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
 polyord = 8;
if nargin > 2 &&  isnumeric(varargin{1})
    polyord = varargin{1};
elseif nargin > 2 
    
    smoothing_method = varargin{1};
    switch smoothing_method
        case 'moving_average'
            if nargin < 4,
                smoothbw = 10;
            else
                smoothbw = varargin{2};
            end
        case 'polynomial'
            if nargin >= 4 && ~isempty( varargin{2})            
                polyord = varargin{2};
            end
    end
end
if nargin > 4
      use_time = varargin{3}; 
else
    use_time = true(size(bx.time));
end

mn = mean(abs(bx.blrep(use_time,:,:)),1);
kp = mn~=0 & usepts;
switch smoothing_method
    case 'polynomial'
        for k = 1:size(mn,3)
            p = polyfit(find(kp(:,:,k))./length(mn(:,:,k)),log(mn(:,kp(:,:,k),k)),polyord);
            bl(:,:,k) = polyval(p,(1:length(mn(:,:,k)))/length(mn(:,:,k)));
        end
    case 'moving_average'
        g = rectwin(smoothbw)';
        bl = log(convn(mn.*kp,g,'same')./convn(kp,g,'same'));
        
end



out = bx.blrep./repmat(exp(bl),size(bx.blrep,1),1);

