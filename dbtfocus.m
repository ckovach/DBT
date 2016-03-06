
function [blout,w0,time,dbif] = dbtfocus(dbx,focus_factor, upsample_factor)

%
% [tf,freq,time,IF] = dbtfocus(dbx,focus_factor,upsample_factor)
%
% DBTFOCUS uses information about instantaneous frequency within
% each band to obtain sub-bandwidth frequency resolution for a DBT
% time-frequency decomposition. Upsampled bands are the original TF 
% amplitude weighted by a gaussian centered at the instantaneous frequency,
% whose width is dbx.bandwidth./focus_factor. 
%
% Input arguments:
%
%     dbx:  A dbt object. Note that the option 'centerDC', must be true.
%
%     focus_factor: Degree of sharpening along the frequency axis.
%
%     upsample_factor: how much to upsample along the frequency axis.
%
%  Output arguments:
%
%     blout: Average band-limited power weighted according to distance from instantaneous frequency. 
%
%     freq: Frequency axis
%
%     time: time axis
%
%     IF: Instantaneous frequencies for each band as time series. 
%
%
% Example:  
%   
%  
%   x = randn(1e4,1) + chirp((1:1e4)'/1e3,10,5,50,'logarithmic') + cos(2*pi*(1:1e4)'/1e3*50);
%   dbx = dbt(x,1000,10,'centerDC',false,'upsampleFx',4,'lowpass',100);
%   [tf,freq,time,IF] = dbtfocus(dbx,4,2);
%   figure, subplot(2,1,1)
%   dbx.specgram;
%   caxis([-15 15]);
%   title('DBT time frequency plot')
%   ylim([0 100])
%   ax = axis;
%   subplot(2,1,2)
%   imagesc(time,freq,20*log10(tf'))
%   caxis([-15 15]);
%   axis xy
%   axis(ax)
%   title('DBTFOCUS output')
%   hold on, 
%   plot(time,IF,'k:')
%   legend('Instantaneous frequency')
% 
% 
% See also DBT and BLPHASE

if nargin < 3
   upsample_factor = 5;   
end

 wres = dbx.bandwidth./(upsample_factor*(1+dbx.upsampleFx));
 
if nargin < 2 || isempty(focus_factor)
    focus_factor = 4;
end

% Use band-limited phase
use_blphase = true;


w = 0:wres:dbx.bandwidth*(1+dbx.shoulder); % demodulated frequencies
w0 = (dbx.bands(1):wres:dbx.bands(end));

win = sin((0:length(w)-1)./length(w)*pi);

%F = ifft(eye(length(dbx.time),length(w)));
arg =@(x)atan2(imag(x),real(x));

if use_blphase
    % Suppress spectral leakage distortion of phase
    dbphn = blphase(dbx);
else
    % No bias corretion
    dbphn = dbx; 
    dbphn.blrep = dbx.blrep./abs(dbx.blrep);
end

blout=zeros([length(dbx.time),length(w0),size(dbx.blrep,3)]);


%%% Get instantaneous frequency.
dbif = arg(dbphn.blrep(1:end,:,:).*conj(dbphn.blrep([1 1:end-1],:,:)))/pi*dbx.sampling_rate/2;

for k = 1:length(dbx.frequency)
        
    wi = find(w0>=dbx.bands(k,1) & w0<dbx.bands(k,2));
    
    
    %%% Gaussian weighting 
    IFmap = exp(-((repmat(w(1:length(wi)),[size(dbif,1),1,size(dbif,3)])-repmat(dbif(:,k,:),1,length(wi))).*focus_factor./dbx.bandwidth).^2);
    
    blout(:,wi,:) = blout(:,wi,:)...
         +  IFmap.*repmat(abs(dbx.blrep(:,k,:)),1,length(wi)).*repmat(win(1:length(wi)).^2,[length(dbx.time) 1 size(dbx.blrep,3)]);
end

time = dbphn.time;

if nargout > 3
    dbif = dbif + repmat(ones(size(dbif,1),1)*dbx.bands(:,1)',[1 1 size(dbif,3)]);
end
