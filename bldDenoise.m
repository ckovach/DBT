function [xdn,F,blsig] = bldDenoise(x,fs,bandwidth,varargin)

% Denoise using the band limited power representation of a signal (see
% BLD).

% C Kovach 2013
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 3 || isempty(bandwidth)
   bandwidth = .05; % This seems to work well for constant line noise. For FM line noise try .25 or so. 
end

use_bld = true;
kurtosis_threshold = 10; % Apply a lower Z threshold to frequency points that have kurtosis above this value

shoulder  = 1; %This is overridden if passed as an argument in varargin

if use_bld
    blsig = bld(x,fs,bandwidth,'padding','time','shoulder',shoulder,varargin{:}); % Band limited representation of the signal (see bld)
else
    blsig = stft(x,fs,bandwidth,'shoulder',shoulder,varargin{:}); % Band limited representation of the signal (see bld)    
end
w = blsig.frequency;

kt = kurtosis(abs(blsig.blrep));
nsig = rmbaseline(blsig,w>=40 & kt<kurtosis_threshold); %takes out the baseline by fitting a polynomial

mn = mean(abs(nsig));


%%% Compute z-score for mean power
nzsc = @(x)(x-nanmean(x))./nanstd(x);
z = nzsc(mn);

%%% exclude points that exceed a threshold from the score computation
z(w<40) = nan;
z(kt>kurtosis_threshold) = nan;
while any(abs(z)>5)
    
    z(abs(z)>5) = nan;
    
    z = nzsc(z);
    
end
    
z = (mn-mean(mn(~isnan(z))))./std(mn(~isnan(z)));
ln = z>5;  %%% Threshold the adjusted z

P = abs(nsig);
P(:,w<40) = nan;
P(:,ln) = nan;

%Compute zscore over all time-frequency points without including the potentially contaminated frequencies in the variance estimate
Z = (abs(nsig)-mean(abs(nsig(~isnan(P)))))./std(abs(nsig(~isnan(P))));

% Exclude frequencies below 40
Z(:,w<40) = 0;
%Set threshold 
LN = isnan(P) & Z > 3 | Z >6 ;

% Smooth edges a little to reduce time-domain artifacts 
g = gausswin(ceil(.5./blsig.sampling_rate));
g = g./max(g);
LN = convn(LN,g./sum(g),'same');

F = 1-LN;

blsig.blrep = blsig.blrep.*F;

xdn = blsig.signal;

xdn = xdn(1:length(x));



