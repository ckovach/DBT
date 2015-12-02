function [xout,spike] = spikefilter(xin,fs,spike)


% Simple script to remove outliers through iterative thresholding. Outliers
% are removed by first applying a z-score transform and rejecting all
% outliers that exceed a given threshold magnitude (defualt is 10), then 
% repeating the first step with remaining data points until no outliers 
% remain beyond the threshold. Data are windowed by a function 
%
%        f = 10^( -A/20*( h(t) o I(t) ) )
%
% where I is the indicator for discarded times, and by default, h(t) is a .2 s
% Hann window and A = 50.
%
% Use: [xfilt,out] = spikefilter(x,fs,[options])
%
% Inputs:
%
%   x - A column vector with data
%   fs - sampling frequency
%   options  - A structure with the following fields: 
%          .threshold  - z score threshold for spikes      
%          .smoothwindow - duration of Hann window if scalar, otherwise the window function
%          .interpolate - option to interpolate:
%                           false: (default) simple windowing, no interpolation. The result output is then
%                                       Xfiltered = X*f 
%                           true: use a Hann window of twice the duration of the smoothing window to interpolate
%                           [scalar value]: use a Hann window of the given duration to interpolate
%                           [vector]: use the given window to interpolate
%               Interpolation computes a moving weighted average, weighted by the window f to give Xinterp, 
%               that is
%                            Xinterp = ( g o (X*f) )/( g o f )  
%               where g is the interpolation window, and the filtered result is  
%                           Xfiltered = X*f + (1-f)*Xinterp

if nargin < 3
    spike.threshold = 10; % This is the threshold used in detecting spikes. 
                          % Z-score is computed iteratively until
                          % no points exceed this value. 
                          
    spike.smoothwindow = .2;% Apply Hann window of given duration to smooth the spike filter.
    
    spike.interpolate = false; % If true, interpolates values through a weighted average, 
                               % with a Hann interpolation window twice
                               % the size of the smoothing window.
                               % If false, then data are subjected to
                               % simple windowing, with excluded points
                               % scaled towards zero.
                               % If a scalar value is given, then uses a
                               % Hann window of the given duration. 
                               % A vector is taken directly as the interpolation window.
                               % If a string, then it uses the matlab
                               % interp1 function with the specified method
                               % (eg. 'spline').
                               
   spike.combine_channels =false; % If true, combine the filter across channels and apply the same filter to all channels.                              
end

DB_attenuation = 50; %Peak attenuation for an isolated impulse

xout = zeros(size(xin));
if size(xin,2)<2
   spike.combine_channels = false; 
end
for k = 1:size(xin,2)
    x = xin(:,k);
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
    
    if isscalar(spike.interpolate) && spike.interpolate
        if islogical(spike.interpolate)
            %%% Use a default interpolation window with twice the support of the smoothing
            %%% window.
            interpwin = hanning(round(length(win)*2));
        else
            %%% Scalar values are treated as the interpolation window duration
            interpwin = hanning(ceil(spike.interpolate.*fs));
        end
    elseif ~isscalar(spike.interpolate)
        %%% For a non-scalar value, the input is the window.
        interpwin = spike.interpolate;
    else
        interpwin = 0;
    end
    
    
    
%     spike.filter(:,k) = exp(convn(log(1-spks+eps),win,'same'));
    spike.filter(:,k) = 10.^(-DB_attenuation/20*convn(spks,win,'same'));
    spike.filter(spike.filter<0,k)=0;
end

if spike.combine_channels
   spike.filter = repmat(prod(spike.filter,2),1,size(xin,2)); 
end

if ~ischar(spike.interpolate)
    %interpwin = win;
    if ~isequal(interpwin,0)
        %%% Smooth x through weighted averaging.
        xconv = convn(xin.*spike.filter,interpwin,'same')./(convn(spike.filter,interpwin,'same')+eps);
        xinterp = + (1-spike.filter).*xconv;
    else
        xinterp = 0;
    end

    xout = xin.*spike.filter + xinterp;
else
    t = (1:size(x,1))';
    xout = interp1(t(spike.filter>.5),xin(spike.filter>.5,:),t,spike.interpolate);
end
    
