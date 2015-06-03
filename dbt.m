classdef dbt

% This class implements the demodulated band transform (DBT). 
% 
% Use:
%
%  B = dbt(X,Fs,BW)
%
%  X  - Signal as a column vector. If X is a matrix, each column is a separate
%        signal.
%  Fs - Sampling frequency
%  BW - Bandwidth of decomposition
% 
%
%  B = dbt(X,Fs,Bw, ['option'], [value])
%
%  B is an object with the following properties:
%
%  B.blrep: time-frequency coefficients with rows corresponding to time
%           and columns to frequency. If X is a matrix, channels correspond to the
%           3rd dimension.
%   .bandwidth:  bandwidth of the DBT.
%   .sampling_rate: sampling rate after the transform.
%   .time:  sampled time points.
%   .frequency: sampled center frequencies/
%   .bands: band limits for each frequency band/
%   .taper: taper used to window frequency bands.
%   .padding: whether signal duration is adjusted through time padding
%             ('time') or fequency padding ('frequency').
%   .fftpad: degree of oversampling achieved through fft padding. 1
%            corresponds to 2x oversampling.
%   .centerDC: If false (default), the fft of the DBT bands contains
%               positive frequencies only, implying 2x oversampling,
%               otherwise each band is demodulated to be centered on DC.
%   
%
%  Options:
%
%  	offset   -  offset of the first band from 0hz (default = 0)
%       padding  - 'time': pad signal in the time domain, changing duration
%                  'frequency': (default) pad in the frequency domain, changing sampling rate
%       shoulder - (0 - 1) degree of overlap between neighboring bands (default = 1). 
%                  Note that, at present, 1 is the maximum allowable value.
%
%  The overlapping portions of the bands are windowed with a  taper:
%	     ____   __________   ________   _________
%	         \ /	      \ /        \ /   
% .  .  .	  X            X          X           . . . 
% 	         / \	      / \        / \
%            |-----------|           |-|
%    	          BW              shoulder
%      
% By default,  the taper is defined so that squares sum to 1
%
% Methods:
%
% B.signal: Returns the signal reconstructed from DBT coefficients via the
%           inverse DBT transform. B.signal(bands) returns the reconstructed 
%           signal including only the specified bands. bands can be a
%           logical or numerical index into frequency bands. 
%           Example:  B.signal(B.frequency < 200) returns the signal including 
%           only bands with center frequency below 200.
%           B.signal(bands,true) returns the complex-valued analytic signal, 
%           containing only positive frequencies.
%           
% B.specgram: plots a spectrogram based on the DBT. B.specgram(true) normalizes
%             the spectrogram with smoothing. 
%
% See also TAPER, DBTCOH, DBTDENOISE, DBTPAC, STFT

%     C Kovach 2013 - 2015
% 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

    properties 
        
        blrep = [];  % Band-limited analytic representaiton of the signal (time x frequency band)
        bandwidth = [];  %The bandwidth parameter. This is the bandwidth not including the shoulders. Full bandwidth is BW*(1+shoulder)
        sampling_rate = []; % Sampling rate after downsampling
        offset = 0;	% Offset if the first band doesn't begin at 0 (a highpass cutoff)
        time = [];     % Time points for the rows in blrep
        frequency = [];  % Center frequencies for the columns in blrep
        bands =[];     % Bands for the columns of blrep 
        fullN = 0;    % Length of the signal after padding
        fullFS = 0;   % Sampling frequency of the original signal
        Norig = 0;    % Original signal length before padding
        shoulder = 1; % Degree of frequency overlap between neighboring bands (0 - 1)
        lowpass = []; % Lowpass cutoff
        taper = [];
    	 fftpad = 0; % Upsample each band by padding its fft by this proportion of the unpadded length.
        userdata=[];
        centerDC = false; % If true, bands are demodulated to [-bandwidth bandwidth], otherwise to [0 bandwidth] with zero padding to 2*bandwidth.
                          % If true, the pass band is centered at 0, so
                          % that the signal contains negative
                          % frequencies, which halves the number of
                          % samples.
                          
         bwtol = 1e-8;    % Tolerance for the bandwidth. Higher values set the bandwidth more precisely but require more padding.           
        direction = 'acausal'; % acausal (default),'causal', or 'anticausal'. Note these are only approximate as strictly causal or anticausal filters can have no zeros on the unit circle.                  
        remodphase = true; % If true, applies a phase correction equivalent to remodulating the subsampled data to the original band. This is necessary to get a
                            % correct average when averaging the spectrum
                            % over time.
          padding = 'time'; % Options are 'time' or 'none'; 'frequency' is obsolete.
    end    
    
    methods
    
        function me = dbt(varargin)
            
            i = 4;       
            me.taper = taper; %#ok<CPROP>
            
           while i < length(varargin)
              switch lower(varargin{i})
                  
                  case {'offset','highpass'}
                      me.offset = varargin{i+1};
                      i = i+1;
                  case 'band'
                      band = varargin{i+1};
                      me.offset = band(1);
                      me.lowpass = band(end);                      
                      i = i+1;
                  case 'padding'
                      me.padding = varargin{i+1};
                      i = i+1;
                  case 'shoulder'
                      me.shoulder = varargin{i+1};
%                       if me.shoulder ~=1
%                           me.shoulder = 1;
%                           warning(sprintf('\n!!!!!!!!!!\nMPORTANT MESSAGE: Due to a bug, shoulder is forced to be 1. This will be fixed in a future revision.\n!!!!!!!!!!'))
%                       end
                      i = i+1;
                  case 'lowpass'
                      me.lowpass = varargin{i+1};
                      i = i+1;
                 case 'taper'
                      me.taper = varargin{i+1};
                      i = i+1;
                 case 'centerdc'
                      me.centerDC = varargin{i+1};
                      i = i+1;
                  case 'fftpad' % Pad fft by proportion of window size
                      me.fftpad = varargin{i+1};
                      i = i+1;
                  case 'upsample' %Upsample by proportion with fft padding
                      me.fftpad = varargin{i+1}-1;
                      i = i+1;
                  case 'direction'
                      me.direction = varargin{i+1};
                      i=i+1;
                  case 'remodphase'
                      me.remodphase = varargin{i+1};
                      i=i+1;
                  case 'bwtol'
                      me.bwtol = varargin{i+1};
                      i=i+1;
                   otherwise
                     error('Unrecognized keyword %s',varargin{i})
              end
              i = i+1;
           end
               
           if isempty(varargin)
               return
           end
           
           fs  = varargin{2};
           bw =  varargin{3};
           
           fullsig  = varargin{1};
           
           n = size(fullsig,1);
           ncol = size(fullsig,2);
           
           me.Norig = n;
           %Resample signal so that everything factors
           %Keeping signal duration fixed and allowing bandwidth and sampling
           %rate to vary.
           
           % T = n./fs; % old signal duration
           %%% K and M need to be integers
            
           
           switch lower(me.padding)
               case {true,'time'}
                   me.padding = 'time';
                   % do nothing.              
               case 'frequency'
                   warning('Frequency padding option is obsolete. Using default time instead')
                   me.padding = 'time';
               case {'none',false}
                   me.padding = 'none';
                   me.bwtol = Inf;
           end
                   
            %%% Pad signal in time so that bandwidth is approximately
            %%% divides padded duration duration
            

            [~,den] = rat(bw/fs/2,me.bwtol);

           
           
           newn = ceil(n/den)*den;
           newT = newn/fs;
           
           winN = round(bw*newT);
           newbw = winN/newT;
           

           me.fullN = newn;
           me.fullFS = fs;
           
           nnyq = ceil((newn+1)/2);
           if isempty(me.lowpass)
              me.lowpass = nnyq/newn*fs ; 
           end
           
           
           %%% Pad the signal
           fullsig(end+1:newn,:) = 0;
           
           F = fft(fullsig./sqrt(newn));
           
           me.bandwidth = newbw;
           
           nsh = min(ceil(me.shoulder*newbw./fs*newn),winN);
           me.shoulder = nsh*fs./newn./newbw;
           
           noffset = round(me.offset*newT -nsh );
           me.offset = noffset*fs/newn + newbw*me.shoulder;
           me.bands(:,1) = (me.offset-newbw*(1+me.shoulder)/2:newbw:me.lowpass-newbw*(1-me.shoulder)/2);%-newfs/newn;
           me.bands(:,2) = me.bands(:,1)+newbw*(1+me.shoulder);%+newfs/newn;
           nwin =size(me.bands,1);

           %%% Reshaping matrix. This step includes the initial
           %%% circular shift.
           rsmat = noffset + repmat((0:nwin-1)*(winN),round(winN*(1+me.shoulder)),1) + repmat((1:round(winN*(1+me.shoulder)))',1,nwin);% -nsh;
           rsmat = mod(rsmat-1,newn)+1;
           dcindx = find(rsmat==1);

             tp = me.taper.make((0:1:nsh-1)/nsh); 
            invtaper = me.taper.make(1-(0:1:nsh-1)/nsh);

           switch me.direction
            %%% Approximate causal or anti-causal filters while
            %%% preserving summation properties of the tapers.
            %%% Note that a truly causal filter cannot have a zero frequency response anywhere.   
               case {'causal','anticausal'}
                   htptp = hilbert([tp,invtaper]);
                   htptp = htptp./(abs(htptp)+eps);
                   if strcmp(me.direction,'causal')
                      htptp = conj(htptp);  
                   end
                   tp = htptp(1:length(tp)).*tp;
                   invtaper = htptp(length(tp)+(1:length(invtaper))).*invtaper;                       
               case 'acausal'
                   % do nothing
               otherwise
                   error('Unrecognized filter direction, %s.',me.direction)
           end

           Frs = zeros([size(rsmat),ncol]);
           for  k = 1:ncol
                f = F(:,k);
               Frs(:,:,k) = double(f(rsmat));
               if nsh>0
                    Frs(end+(1-nsh:0),1:nwin,k) = diag(sparse(tp))*Frs(end+(1-nsh:0),1:nwin,k);

                   Frs(1:nsh,1:nwin,k) = diag(sparse(invtaper))*Frs(1:nsh,1:nwin,k); 
               end

                %%% Set all negative frequencies to zero
%                     Frs(rsmat(:, end-1:end)>ceil((newn+1)/2),end,k)=0;
%                     Frs(rsmat(:, 1:2)>ceil((newn+1)/2),1,k)=0;

                %%% Corrections for DC and Nyquist to preserve power
                %%% after zeroing negative frequencies
%                 Frs(rsmat(:,1)==1,1,k)=Frs(rsmat(:,1)==1,1,k)/sqrt(2); 
                Frs(dcindx + (k-1)*numel(rsmat))=Frs(dcindx + (k-1)*numel(rsmat))/sqrt(2); 

           end

           padN = floor((me.fftpad+1)*winN*(1+~me.centerDC)*(1+me.shoulder));
           me.fftpad = 1/(1+me.shoulder)*padN/winN/(1+~me.centerDC)-1;
           
           if padN > round(winN*(1+me.shoulder))
               Frs(padN,:,:) = 0;
           elseif padN < round((1+me.shoulder)*winN)
               padN =round( (1+me.shoulder)*winN);
               me.fftpad = 0;
           end
           
           if me.centerDC
               Frs = circshift(Frs,-ceil(winN/2*(1+me.shoulder)));
               me.blrep = ifft(Frs)*sqrt(padN)*sqrt(2);
           else

               me.blrep = ifft(Frs)*sqrt(padN)*sqrt(2);
           end
           
               
            me.sampling_rate = padN/newT;
                 
         
           me.time = (((1:size(me.blrep,1))-1)*newT./size(me.blrep,1))';
           w = mean(me.bands,2)';
           me.blrep(:,w > me.lowpass + me.shoulder*bw,:) = [];
           w(w > me.lowpass + me.shoulder*bw) = [];
           me.frequency = w;

           if me.remodphase  
               %%% This remodulates the signal to obtain the correct phase
               %%% at a given frequency for the original band. It is necessary 
               %%% when averaging the complex-valued spectrum over time so 
               %%% that phase aligns correctly.
               remodulator = exp(1i*2*pi*me.time*(me.frequency - (~me.centerDC)*me.bandwidth*(1+me.shoulder)/2));
               me.blrep = me.blrep.*repmat(remodulator,[1 1 size(me.blrep,3)]);
               
           end
           
        end

        %%%%

        function [data,fs] = signal(me,columnfilter,doHilbert,return_padded)
            
            % db.signal
            % Reconstruct the signal from its band-limited analytic
            % representation by applying the inverse DBT. 
            %
            % db.signal(bands)
            %
            % Reconstruct the signal using only the specified bands. 'bands
            % may be a numeric or logical index. 
            % Example: 
            %        db.signal(db.frequency > 100)
            % returns the signalusing only bands with center frequency above 100 Hz.
            %
            % db.signal(bands,true)
            %
            % Returns a complex-valued analytic signal with negative frequencies zero'd. 
            %
            
            
            
            if nargin < 2 || isempty(columnfilter)
                mult = 1;
            elseif islogical(columnfilter)
                mult = diag(sparse(columnfilter));
            elseif min(size(columnfilter)) == 1 && min(columnfilter)>=1
                
               mult = diag(sparse( ismember(1:size(me.blrep,2),columnfilter)));
               
            elseif min(size(columnfilter)) == 1
                mult = diag(sparse(columnfilter));
            else
                mult = columnfilter;
            end
            
            if nargin < 3 || isempty(doHilbert)
               doHilbert = false; 
            end
            if nargin < 4 || isempty(return_padded)
               return_padded = false; 
            end
            n = me.fullN;
            noffset = round((me.offset-me.bandwidth*me.shoulder)./me.fullFS*n);
            ncol = size(me.blrep,3);
           if me.remodphase  
               %%% If phase remodulation was applied we need to reverse it.
               demodulator = exp(-1i*2*pi*me.time*(me.frequency - (~me.centerDC)*me.bandwidth*(1+me.shoulder)/2));
               me.blrep = me.blrep.*repmat(demodulator,[1 1 size(me.blrep,3)]);

           end

            F = zeros(size(me.blrep));
            for k = 1:ncol
                F(:,:,k) = fft(me.blrep(:,:,k) )*mult/sqrt(size(me.blrep,1)/2);         
            end
            nsh = round(me.shoulder*me.bandwidth./me.fullFS*me.fullN);
            padN = length(me.time);
            winN = round(padN./(1+me.fftpad)/(1+~me.centerDC)/(1+me.shoulder)); 
            if me.centerDC
               
                F = circshift(F,ceil(winN/2*(1+me.shoulder)));%sqrt(2);
                
            end
            F = F(1:round(winN*(1+me.shoulder)),:,:)*sqrt(2);

             nnyq = size(F,1);
            if nsh >0
                % Taper is normally defined so that h(k).^2 + h(k+bw).^2 = 1
                tp = me.taper.make((0:1:nsh-1)/nsh); 
                invtaper = me.taper.make(1-(0:1:nsh-1)/nsh);
                switch me.direction
                %%% Approximate causal or anti-causal filters while
                %%% preserving summation properties of the tapers.
                   
                   case {'causal','anticausal'}
                       htptp = hilbert([tp,invtaper]);
                       htptp = htptp./(abs(htptp)+eps);
                       if strcmp(me.direction,'anticausal')
                          htptp = conj(htptp);  
                       end
                       tp = htptp(1:length(tp)).*tp;
                       invtaper = htptp(length(tp)+(1:length(invtaper))).*invtaper;                       
                    case 'acausal'
                        % Do nothing.
                   otherwise
                       error('Unrecognized filter direction, %s.',me.direction)
                end
               
                for k = 1:ncol
                    sh = diag(sparse(invtaper))*F(1:nsh,2:end,k);
                    
                    F(nnyq-nsh+1:nnyq,1:end-1,k) = diag(sparse(tp))*F(nnyq-nsh+1:nnyq,1:end-1,k)+sh;
                end
                 
            else
                invtaper = 1;
            end
            
            if nsh > 0
                %%% Retain the leading edge if there is an offset
               F0 = permute(F(1:nsh,1,:),[1 3 2]).*repmat(invtaper',1,ncol); 
            else
                F0=zeros(0,ncol);
            end
            
            F(1:nsh,:,:)=[];
            Ffull = zeros(me.fullN,ncol);
            switch me.padding
                case {'frequency','fft'}
                    warning('Frequency padding is an obsolete option.')
            end

            for k = 1:ncol
                f = F(:,:,k);
  
                Ffull( 1 + mod(noffset+(1:numel(f)+size(F0,1))-1,me.fullN),k) = [F0(:,k);f(:)] * sqrt(me.fullN);

               if ~mod(me.fullN,2)
                  Ffull(ceil(me.fullN/2+1),k) = Ffull(ceil(me.fullN/2+1),k)/2;
               end
               Ffull(1,k) = real(Ffull(1,k))/sqrt(2);
            end
            Ffull(ceil((me.fullN+1)/2)+1:end,:,:) = 0;
            
            Ffull = Ffull./sqrt(2);
            if doHilbert %%% If second argument is true return analytic signal
                data = ifft(full(Ffull));
            else
                data = real(ifft(full(Ffull)));
            end
            if ~return_padded
                data = data(1:me.Norig,:);
            end
            fs = me.fullFS;
            
        end
       %%%%
        function varargout = specgram(me,normalize)
            if nargin < 2
                normalize = 0;
            end
            if normalize
                fun = @rmbaseline;
            else
                fun = @(x)x.blrep;
            end
            S = 20*log10(abs(fun(me)))';
            t = me.time;
            w = me.frequency;
            if nargout == 0
               imagesc(t,w,S) 
               axis xy
            else
                varargout(1:3) = {S,t,w};
            end
        end
    end
end
    
