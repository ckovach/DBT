
function [out,bl] = rmbaseline(bx,usepts,varargin)


% [out,bl] = rmbaseline(bx,usepts,varargin)
%
% Normalizes a time-frequency plot with an estimated baseline spectrum.
%
% Keyword options:
%       'smoothing method' : Method used to estimate the baseline (noise-free)
%                                spectrum (see RMBASELINE) 
%                             'moving average' - Baseline is computed iteratively 
%                                as a local average in the frequency domain, 
%                                excluding points that deviate more than 3 
%                                stdev at each step (Default as of rev. 801 2/27/2017).
%                             'local' - as above except applied within a
%                                local moving time window.
%                              'polynomial' - estimate baseline sepctrum as a 10th
%                                order polynomial (Default before rev. 801).
% 
%       'smoothbw'          : Bandwidth of the frequency smoothing winow
%                             used to obtain the local average for the 'moving
%                             average' and 'local' methods (Default = 10 Hz).
% 
%       'smoothbwn'         : Same as above except expressed as a number of
%                             frequency samples. 
% 
%       'local_time_windown': Number of time samples within the moving
%                             window used for the 'local' option.
% 
%       'smoothing polyord' :  Order of the polynomial used for the
%                               'polynomial' smoothing option (Default = 10).
%
% C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

smoothing_method = 'local';
 polyord = 8;
%%% the following apply to moving average and local smoothing
smoothbwn=0;
smoothbw = 10;
smoothwin = 'hann';
smooth_threshold=3;%Iterative discard points with z-score above this threshold from smoothing.
%%% the following apply to local
local_time_smoothingN = 5; %apply a second smoothing over this many time samples in setting the baseline threshold.
local_time_smoothing_win = 'hann';

use_time = true(size(bx.time));
i = 1;

while i <= length(varargin)
    switch lower(varargin{i})
        
       case {'smoothing method'}  
            smoothing_method = varargin{i+1};                       
            i = i+1;
        case 'polyord'
            polyord = varargin{i+1};                    
            i = i+1;
         case {'smoothing bandwidth','smoothbw'}
             smoothbw = varargin{i+1};                    
             i = i+1; 
        case 'smoothing threshold'
            smooth_threshold = varargin{i+1};                    
            i = i+1; 
        case 'local_time_windown'
            local_time_smoothingN = varargin{i+1};                    
            i = i+1; 
        case {'smooth bwn','smoothing bwn'}
            smoothbwn = varargin{i+1};                    
            i = i+1; 
        case 'use time'
            use_time = varargin{i+1};                    
            i = i+1; 
       otherwise
            if ischar(varargin{i})
                error('Unrecognized keyword %s',varargin{i});
            else
                error('Keyword must be a string')
            end
    end
    i = i+1;
end        

if nargin < 2 || isempty(usepts)
    usepts = true;
end
% if nargin > 2 &&  isnumeric(varargin{1})
%     polyord = varargin{1};
% elseif nargin > 2 
%     
% %     smoothing_method = varargin{1};
%     switch smoothing_method
%         case 'moving_average'
%             if nargin < 4,
%                 smoothbw = 10;
%             else
%                 smoothbw = varargin{2};
%             end
%         case 'polynomial'
%             if nargin >= 4 && ~isempty( varargin{2})            
%                 polyord = varargin{2};
%             end
%     end
% end
% if nargin > 4
%       use_time = varargin{3}; 
% else
%     use_time = true(size(bx.time));
% end
B = abs(bx.blrep(use_time,:,:));
switch lower(smoothing_method)
    
    case 'polynomial'
        B = B + 0./B;
        mn = exp(nanmean(log(B),1));
        kp = mn~=0 & usepts;    p = polyfit(find(kp)./length(mn),log(mn(kp)),polyord);
        bl = polyval(p,(1:length(mn))/length(mn));        
        BL = repmat(exp(bl),size(bx.blrep,1),1);
    case {'moving_average','moving average','local'}
        if smoothbwn ==0
            smoothbwn = ceil(smoothbw./bx.bandwidth*(bx.upsampleFx+1));
        end
        g = window(smoothwin,smoothbwn)';
        g = g./sum(g);
        switch lower(smoothing_method)
            case {'moving_average','moving average'}
                B = B + 0./B;
                mn = exp(nanmean(log(B),1));
                kp = mn~=0 & usepts;

                bl = convn(log(mn).*kp,g,'same')./convn(kp,g,'same');
                if smooth_threshold ~=0 && ~isinf(smooth_threshold)

                     thr = iterz(bl-log(mn),smooth_threshold);
                     kp = kp & ~isnan(thr);
                     bl = convn(log(mn).*kp,g,'same')./convn(kp,g,'same');
                end                
                BL = repmat(exp(bl),size(bx.blrep,1),1);
            case 'local'
                B(isnan(B))=0;
                kp = B~=0;
                if length(usepts)>1
                    kp = kp& repmat(usepts,size(B,1),1);
                end
                twin = window(local_time_smoothing_win,local_time_smoothingN);
                twin = twin./sum(twin);
                MN = convn(log(B),twin(:),'same');
                KP = convn(kp,g,'same');
                
                bl = convn(MN.*KP,g,'same')./convn(KP,g,'same');
                if smooth_threshold ~=0 && ~isinf(smooth_threshold)

                     thr = iterz(MN'-bl',smooth_threshold)';
                     kp = kp & ~isnan(thr);
                     bl = convn(MN.*kp,g,'same')./convn(kp,g,'same');
                end                
                BL = zeros(size(bx.blrep));
                BL(use_time,:) = exp(bl);
        end
    otherwise
        error('Unrecognized smoothing method.')
       
end


out = bx.blrep./BL;

