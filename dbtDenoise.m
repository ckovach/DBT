function [xdn,F,blsig,spike] = dbtDenoise(x,fs,bandwidth,varargin)

% Denoise using the demodulated band representation of a signal (see
% DBT). This works by applying a threshold to coefficients of the 
% decomposition, discarding those whose magnitude exceeds the threshold
% and reconstructing the denoised signal from the inverse DBT transform. 
% The threshold is adjusted according to the location of sharp peaks in the
% periodogram over the whole recording and the kurtosis within each band.
% 
%
% Usage:
%
%   [xdn,filt,dbx] = dbtDenoise(x,Fs,bandwidth,[keyword],[value])
%
%    Inputs:
% 
%       x - signal
%       Fs - sampling frequency
%       bandwidth to use in the dbt (default = 0.05)
%
% 
%    Outputs:
% 
%       xdn   - denoised signal
%       filt  - filter used on DBT coefficients in denoising  
%       dbx - DBT of xdn
%
%    Keyword options:
%       
%       'remove spikes':  true - remove spikes at the default threshold (default)
%                                before applying dbt-based denoising.
%                         false - no spike removal
%                         n (scalar) - remove spikes at specified threshold
%                             Thresholding is applied iteratively to z-scores
%                             across all remaining time points none remain above 
%                             the threshold (see SPIKEFILTER).
%
%       'kthresh'      :  kurtosis threshold. Coefficients in frequency bands for which
%                             kurtosis over time falls above this value are
%                             thresholded at 'zlothreshold' and otherwise at
%                             'zhireshold'. This feature was added as
%                             kurtosis is useful for detecting frequency
%                             modulated line noise. Default = 10.
%                         
%      'zhithresh'     : Threshold for coefficient removal. (default=6)
%       
%      'zlothresh'     : Threshold for coefficient removal in frequency
%                        banbds that fall above 'kthresh'. (default = 3).
%       
% See also DBT, SPIKEFILTER

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2013

if nargin < 3 || isempty(bandwidth)
   bandwidth = .05; % This seems to work well for constant line noise. For FM line noise try .25 or so. 
end

kurtosis_threshold = 10; % Apply a lower Z threshold to frequency points that have kurtosis above this value
spike.remove_spikes = true;  % Zero out time points that exceed some threshold
spike.threshold = 10; % This is the threshold used in detecting spikes. 
                      % Z-score is computed iteratively until
                      % no points exceed this value. The threshold is set
                      % high by default because the main purpose here is to avoid
                      % distortion of kurtosis used in the
                      % kurtosis-threshold.
spike.smoothwindow = .2;% Apply hanning window of given duration to smooth the spike filter.
spike.interpolate = false;% 
filter_above = 40; 
use_stft = false;
zhithresh = 6;
zlothresh = 3;
flagthresh = 3;
makeplots = false;
smoothing_method = 'polynomial';
adjust_threshold = true; % If true this performs an initial denoising run at a higher threshold if an excessive number of frequency bands are rejected (>15 %)
prefilter_threshold = 10; % Percent rejected bands needed to trigger prefiltering at a higher threshold
baseline_polyord = 10;
smbw = 10;
rm_edge_samples = 2; %Remove this many edge samples 
i = 1;

argin = varargin;

while i <= length(varargin)
          switch lower(varargin{i})

           case {'dbt'}  % use stft instead of dbt
                  use_stft = false;                        
                  varargin(i) = [];
                  %i = i-1;
              case {'kurtosis','kthresh'}  % Apply a lower filter threshold for frequency values that have kurtosis exceeding this value
                  kurtosis_threshold = varargin{i+1};
                  varargin(i:i+1) = [];
                  i = i-1;
               case {'filter above'}  % ONly apply filter above this frequency
                  
                  filter_above = varargin{i+1};
                  varargin(i:i+1) = [];
                  i = i-1;
               case {'remove spikes'}  
                  
                  spike.remove_spikes = varargin{i+1};
                  if ~islogical(spike.remove_spikes) && spike.remove_spikes~=1 && spike.remove_spikes~=0
                     spike.threshold = spike.remove_spikes; 
                  end
                  varargin(i:i+1) = [];
                  i = i-1;
                 case {'spike window'}  % spike window widht 
                    spike.smoothwindow = varargin{i+1};
                    varargin(i:i+1) = [];
                     i = i-1;
                 case {'spike opts'}  % option structure for spike exclusion (see SPIKEFILTER)
                    spike = varargin{i+1};
                    varargin(i:i+1) = [];
                    if ~isfield(spike,'remove_spikes')
                        spike.remove_spikes = true;
                    end
                     i = i-1;
                case {'zhithresh','high threshold'}  % coefficient threshold
                  zhithresh = varargin{i+1};
                  varargin(i:i+1) = [];
                  i = i-1;        
                case {'zlothresh','low threshold'}  % Apply a lower threshold to frequencies above the kurtosis threshold
                  zlothresh = varargin{i+1};
                  varargin(i:i+1) = [];
                  i = i-1;   
                case {'flagthresh','flag threshold'}  % Threshold baseline-corrected z score at which to identify bands as contaminated.
                  flagthresh = varargin{i+1};
                  varargin(i:i+1) = [];
                  i = i-1;   
                case {'makeplots'} 
                  makeplots = varargin{i+1};
                  varargin(i:i+1) = []; 
                  i = i-1;                
                case {'smoothing method'} 
                  smoothing_method = varargin{i+1};
                  varargin(i:i+1) = []; 
                  i = i-1;               
  
                case {'smoothing bandwidth'} 
                  %Bandwidth of moving average for 'moving average' method
                  smbw = varargin{i+1};
                  varargin(i:i+1) = []; 
                  i = i-1;               
                case {'smoothing polyord'} 
                  baseline_polyord = varargin{i+1};
                  varargin(i:i+1) = []; 
                  i = i-1;               
                case {'adjust threshold'} 
                  adjust_threshold = varargin{i+1};
                  if ~islogical(adjust_threshold)
                      prefilter_threshold = adjust_threshold;
                      adjust_threshold = true;
                  end
                  
                  varargin(i:i+1) = []; 
                  i = i-1; 
             case {'remove edge','rm edge'} 
                  %Number of edge samples to remove
                  rm_edge_samples = double(varargin{i+1});
                  varargin(i:i+1) = []; 
                  i = i-1;  
              otherwise
%                  error('Unrecognized keyword %s',varargin{i})
          end
          i = i+1;
end
  

shoulder  = 1; %This is overridden if passed as an argument in varargin

if spike.remove_spikes    
    [x,spike] = spikefilter(x,fs,spike);
end

if ~use_stft
    blsig = dbt(x,fs,bandwidth,'padding','time','shoulder',shoulder,'upsample',1,varargin{:}); % Band limited representation of the signal (see dbt)
    nsmbw = ceil(smbw./blsig.bandwidth);
else
    blsig = stft(x,fs,bandwidth,'shoulder',shoulder,varargin{:}); % Band limited representation of the signal (see dbt)    
    nsmbw = ceil(smbw.*blsig.timewindow);
end
w = blsig.frequency;

include_times = blsig.time >= rm_edge_samples./blsig.sampling_rate  ...
               & blsig.time <= blsig.Norig./blsig.fullFS - rm_edge_samples./blsig.sampling_rate;

kt = kurtosis(abs(blsig.blrep(include_times,:,:)));

pcntrej = mean(kt>kurtosis_threshold);
F0=1;
if adjust_threshold && pcntrej > prefilter_threshold/100;
    fprintf('\n%0.0f%% bands flagged. Prefiltering with higher rejection threshold.\n    ',pcntrej*100)
    [x,F0,blsig] = dbtDenoise(x,fs,bandwidth,argin{:},'low threshold',2*zlothresh,'flag threshold',2*flagthresh,'adjust threshold',true,'kthresh',2*kurtosis_threshold,'spike opts',spike);
    kt = kurtosis(abs(blsig.blrep(include_times,:,:)));
end


switch smoothing_method
    case 'polynomial'
        nsig = rmbaseline(blsig,w>=filter_above & kt<kurtosis_threshold,smoothing_method,baseline_polyord,include_times); %takes out the baseline by fitting a polynomial
    case 'moving_average'
        nsig = rmbaseline(blsig,w>=filter_above & kt<kurtosis_threshold,smoothing_method,nsmbw,include_times); %takes out the baseline by fitting a polynomial
        
end
mn = mean(abs(nsig(include_times,:,:)));


%%% Compute z-score for mean power
nzsc = @(x)(x-nanmean(x))./nanstd(x);
z = nzsc(mn);

%%% exclude points that exceed a threshold from the score computation
z(w<filter_above) = nan;
z(kt>kurtosis_threshold) = nan;

while any(abs(z)>flagthresh)
    
    z(abs(z)>flagthresh) = nan;
    
    z = nzsc(z);
    
end
    
z = (mn-mean(mn(~isnan(z))))./std(mn(~isnan(z)));
ln = z>flagthresh;  %%% Threshold the adjusted z

P = abs(nsig);
P(:,w<filter_above) = nan;
P(:,ln) = nan;
P(~include_times,:) = nan;
%Compute zscore over all time-frequency points without including the potentially contaminated frequencies in the variance estimate
Z = (abs(nsig)-mean(abs(nsig(~isnan(P)))))./std(abs(nsig(~isnan(P))));

% Exclude frequencies below 40
Z(:,w<filter_above) = 0;
%Set threshold 
LN = isnan(P) & Z > zlothresh | Z >zhithresh ;



F = (1-LN).*F0;


%%% Apply the filter
blsig.blrep = blsig.blrep.*F;

xdn = blsig.signal;

xdn = xdn(1:length(x));


if makeplots
%    dnnsig = blsig.blrep; 
   pl = plot(blsig.frequency,100*(1-mean(F))');
   if ~isequal(F0,1)
        hold on
        pl(2) = plot(blsig.frequency,100*(1-mean(F0))');
        set(pl(2),'color','b');
        hold off
        
        legend({'Final','First pass'})
   end
       
   set(pl(1),'color','r');
   xlabel('Freq (hz)')
   ylabel('% Discarded');
   grid on
   xlim([blsig.frequency(1) blsig.frequency(end)])
end


