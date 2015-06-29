function [xfilt,spike] = spikefilter(x,fs,spike)


% Simple script to remove outliers through iterative thresholding. Outliers
% are removed by first applying a z-score transform and rejecting all
% outliers that exceed a given threshold magnitude (defualt is 10), then 
% repeating the first step with remaining data points until no outliers 
% remain beyond the threshold. Data are windowed by a function that convolves
% a tapered log(1 - Hann + eps) window with the  
%
% Use: [xfilt,out] = spikefilter(x,fs,[options])
%
% Inputs:
%
%   x - A column vector with data
%   fs - sampling frequency
%   options  - A structure with the following fields: 
%          .threshold  -       
%
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 3
    spike.threshold = 10; % This is the threshold used in detecting spikes. 
                          % Z-score is computed iteratively until
                          % no points exceed this value. The threshold is set
                          % high by default because the main purpose here is to avoid
                          % distortion of kurtosis used in the
                          % kurtosis-threshold.
    spike.smoothwindow = .2;% Apply hanning window of given duration to smooth the spike filter.

end

spks = false(size(x));
newspks = true;
while any(newspks)
       z = (x-mean(x(~spks)))/std(x(~spks));
       newspks = abs(z)>spike.threshold &~spks; 
       spks = newspks | spks;
end
if isscalar(spike.smoothwindow)
    win = hanning(ceil(spike.smoothwindow.*fs));
else
    win = spike.smoothwindow; 
end
spike.filter = exp(convn(log(1-spks+eps),win,'same'));
spike.filter(spike.filter<0)=0;
xfilt = x.*spike.filter;