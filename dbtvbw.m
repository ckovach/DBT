classdef dbtvbw

% This class implements the demodulated band transform (DBT) with variable 
% bandwidths. 
% 
% Use:
%
%  B = dbt(X,Fs,BW)
%
%  X  - Signal matrix with channels in columns. 
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
%   .frequency: sampled center frequencies.
%   .bands: band limits for each frequency band.
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
%   padding  - 'time' (default) or 'none': pad signal in the time domain to
%                adjust bandwidth.
%   shoulder - (0 - 1) degree of overlap between neighboring bands (default = 1). 
%                  Note that, at present, 1 is the maximum allowable value.
%
%      The overlapping portions of the bands are windowed with a  taper:
%              ____   __________   ________   _________
%                  \ /	        \ /        \ /   
%     .  .  .	  X              X          X           . . . 
%                  / \	        / \        / \
%                  |-----------|           |-|
%                   BW              shoulder
%      
%       By default,  the taper is defined so that squares sum to 1
%
%  upsampleFx: upsample the frequency scale by a factor of (1+x) (default x
%              = 0, no upsampling ).
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
        upsampleFx = 0; %Amount by which to oversample the frequency dimension.
        userdata=[];
        cfreq = [];
        centerDC = true; % If true, bands are demodulated to [-bandwidth bandwidth], otherwise to [0 bandwidth] with zero padding to 2*bandwidth.
                          % If true, the pass band is centered at 0, so
                          % that the signal contains negative
                          % frequencies, which halves the number of
                          % samples.
        logspace = false;                  
         bwtol = 1e-8;    % Tolerance for the bandwidth. Higher values set the bandwidth more precisely but require more padding.           
        direction = 'acausal'; % acausal (default),'causal', or 'anticausal'. Note these are only approximate as strictly causal or anticausal filters can have no zeros on the unit circle.                  
        remodphase = false; % If true, applies a phase correction equivalent to remodulating the subsampled data to the original band. This is necessary for example to get a
                            % correct average when averaging the raw spectrum over time.
          padding = 'time'; % Options are 'time' or 'none'; 'frequency' is obsolete.
         inputargs = {};
 
    end    
    
    methods
    
        function me = dbtvbw(varargin)
            
            i = 3;       
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
                  case {'cfreq','center frequency','center_frequency'}
                      me.cfreq = varargin{i+1};
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
                  case {'fftpad','upsampletx'} % Pad fft by proportion of window size
                      me.fftpad = varargin{i+1};
                      i = i+1;
                  case 'upsample' %Upsample time by proportion with fft padding and frequency by increasing overlap
                      nx =  varargin{i+1}-1;
                      me.fftpad = nx;
                      me.upsampleFx = nx;
                      if nx ~= round(nx)
                          warning('%s will not currently allow signals to be reconstructed with non-integer upsampling. Integer upsampling is recommended.',upper(mfilename))
                      end
                      i = i+1;
                  case 'upsamplefx' %Upsample frequency by proportion + 1 
                      
                      nx =  varargin{i+1};
                    
                      if nx ~= round(nx)
                        warning('%s will not currently allow signals to be reconstructed with non-integer upsampling. Integer upsampling is recommended.',upper(mfilename))
                      end
            
                      me.upsampleFx = nx;
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
                   case 'logspace'
                      me.logspace = varargin{i+1};
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
%            bw =  varargin{3};
           
           fullsig  = varargin{1};
           
           n = size(fullsig,1);
           ncol = size(fullsig,2);
           
           me.inputargs = varargin(3:end);
        
           me.Norig = n;
           
            nnyq = ceil((n+1)/2);
           if isempty(me.lowpass)
              me.lowpass = nnyq/n*fs ; 
           end
           
           %%% A vector giving the "proportion" of the distance from one center frequency to the next 
           me.cfreq(me.cfreq<me.offset | me.cfreq>me.lowpass)=[];
           
           cfreq0 = unique([me.offset-(2-1./(me.upsampleFx+1)).*(min(me.cfreq)-me.offset); me.offset;me.cfreq(:);me.lowpass; me.lowpass+(me.lowpass-max(me.cfreq)) ]*2/fs);
           
           
           w = mod(2*(0:n-1)./n + 1,2)-1;
           % Create a margin at nyquist for the sake of consistency
           w(w<max(cfreq0)-2) =  w(w<max(cfreq0)-2)+2 ;
       
           %%% Circular shift to simplify things
           ncs = floor(min(cfreq0)*n/2);
           wcs = circshift(w,[0,-ncs]);
       
           for k = 1:me.upsampleFx+2
         
               cfrq = cfreq0(1:end-1) + diff(cfreq0)*(k-1)./(me.upsampleFx+1);            
               frqscale(:,k) = mod(interp1(cfrq,1:length(cfrq),wcs,'linear'),1);
               frqscale(:,k) = mod(interp1(cfrq,1:length(cfrq),wcs,'linear'),1);  
                if me.logspace 
                   cfrq(cfreq0(1:end-1)>0) = exp(log(cfreq0(cfreq0(1:end-1)>0)) + diff(log(cfreq0(cfreq0>0)))*(k-1)./(me.upsampleFx+1));            
                   frqscale(wcs>0,k) = mod(interp1(log(cfrq(cfrq>0)),1:sum(cfrq>0),log(wcs(wcs>0)),'linear'),1);  
                       
                end
               cfreqs(k,:) = cfrq;
           end
           intval = [zeros(1,size(frqscale,2));cumsum(diff(frqscale)<0)];
           tp = me.taper;
           
           Wgt = tp.make(frqscale);
           
           F = circshift(fft(fullsig),[-ncs,0]);           
           
           
           %Apply weighting
           srate = [];
           for k = 1:size(cfreq0)-2
              for kk = 1:size(Wgt,2)
                   
                   indx1 = intval(:,kk)==k-1 & ~isnan(Wgt(:,kk));
                   indx2 = intval(:,kk)==k & ~isnan(Wgt(:,kk));
                   
                   win1 = sqrt(1-Wgt(indx1,kk).^2);
                   win2 = Wgt(indx2,kk);
                   n1 = length(win1);
                   n2 = length(win2);
                   mxn = ceil(max(n1,n2)*(1+me.fftpad));
                   
                   f1 = zeros(mxn,ncol); 
                   f2=f1;
                   f1(end-n1+1:end,:) = F(indx1,:).*repmat(win1,1,ncol);
                   f2(1:n2,:) = F(indx2 ,:).*repmat(win2,1,ncol);
                   f = cat(1,f2,f1);
                   
                   blsigs{kk,k} = ifft(f);
                   times{kk,k} = (0:length(f)-1)./length(f)*n./fs;
                   ws = wcs(indx1|indx2);
                   frange{kk,k} = [min(ws) max(ws)]*fs/2;
                   srate(kk,k) = length(f)./n*fs;
              end
              
           end
          
           me.blrep = blsigs;
           me.time = times;
           me.frequency = cfreqs(:,2:end)*fs/2;
           me.sampling_rate = srate;
           me.bands = frange;
           me.fullFS = fs;
           
        end
        function sig = signal(me,varargin)
            
            
                 
           fs  = me.fullFS;
%            bw =  varargin{3};
           
           
           n = me.Norig;
           ncol = size(me.blrep{1},2);
           
          
           
           cfreq0 = unique([me.offset-(min(me.cfreq)-me.offset); me.offset;me.cfreq(:);me.lowpass; me.lowpass+(me.lowpass-max(me.cfreq)) ]*2/fs);
           
           
           w = mod(2*(0:n-1)./n + 1,2)-1;
           % Create a margin at nyquist for the sake of consistency
           w(w<max(cfreq0)-2) =  w(w<max(cfreq0)-2)+2 ;
       
           %%% Circular shift to simplify things
           ncs = floor(min(cfreq0)*n/2);
           wcs = circshift(w,[0,-ncs]);
       
           for k = 1:me.upsampleFx+2 
               cfrq = cfreq0(1:end-1) + diff(cfreq0)*(k-1)./(me.upsampleFx+1);            
               frqscale(:,k) = mod(interp1(cfrq,1:length(cfrq),wcs,'linear'),1); %#ok<AGROW>
                if me.logspace 
                   cfrq(cfreq0(1:end-1)>0) = exp(log(cfreq0(cfreq0(1:end-1)>0)) + diff(log(cfreq0(cfreq0>0)))*(k-1)./(me.upsampleFx+1));            
                   frqscale(wcs>0,k) = mod(interp1(log(cfrq(cfrq>0)),1:sum(cfrq>0),log(wcs(wcs>0)),'linear'),1);   %#ok<AGROW>
                       
                end
           end
           
           intval = [zeros(1,size(frqscale,2));cumsum(diff(frqscale)<0)];
           tp = me.taper;
           
           Wgt = tp.make(frqscale);

           %Reverse
           F = zeros(me.Norig,ncol);
           for k = 1:size(cfreq0)-2
                for kk = 1:size(Wgt,2)
                   
                   indx1 = intval(:,kk)==k-1 & ~isnan(Wgt(:,kk));
                   indx2 = intval(:,kk)==k & ~isnan(Wgt(:,kk));
                   
                   win1 = sqrt(1-Wgt(indx1,kk).^2);
                   win2 = Wgt(indx2,kk);
                   n1 = length(win1);
                   n2 = length(win2);
                   
                   f = fft(me.blrep{kk,k}); 
                   f1 = f(end-n1+1:end,:);
                   f2 = f(1:n2);
                   F(indx1,:) = f1.*repmat(win1,1,ncol)+F(indx1,:);
                   F(indx2,:) = f2.*repmat(win2,1,ncol)+F(indx2,:);

                end
              
           end
           F(wcs>1 | wcs< 0,:) = 0;
           sig = real(ifft(circshift(F,[ncs,0])));           
        end
        
        function specgram(me)
        
            hold on
            mnf = cellfun(@min,me.bands);
            mxf = cellfun(@max,me.bands);      
            plotmn = cellfun(@(x) x + [1 -1]*diff(x)/3,me.bands);
          
        end
        
%            %Resample signal so that everything factors
%            %Keeping signal duration fixed and allowing bandwidth and sampling
%            %rate to vary.
%            
%            % T = n./fs; % old signal duration
%            %%% K and M need to be integers
%             
%            
%            switch lower(me.padding)
%                case {true,'time'}
%                    me.padding = 'time';
%                    % do nothing.              
%                case 'frequency'
%                    warning('Frequency padding option is obsolete. Using default time instead')
%                    me.padding = 'time';
%                case {'none',false}
%                    me.padding = 'none';
%                    me.bwtol = Inf;
%            end
%                    
%             %%% Pad signal in time so that bandwidth is approximately
%             %%% divides padded duration duration
%             
% %             [~,den] = rat(bw/fs/2,me.bwtol);
%             stepsize = bw/(me.upsampleFx +1);           
%             [~,den] = rat(stepsize/fs/2,me.bwtol);
% 
%            
%            
%            newn = ceil(n/den)*den;
%            newT = newn/fs;
%            
%            winN = round(bw*newT);
%            newbw = winN/newT;
%            
% 
%            me.fullN = newn;
%            me.fullFS = fs;
%            
%           
%            
%            
%            %%% Pad the signal
%            fullsig(end+1:newn,:) = 0;
%            
%            F = fft(fullsig./sqrt(newn));
%            
%            me.bandwidth = newbw;
%            
%            
%            nsh = min(ceil(me.shoulder*newbw./fs*newn),winN);
%            me.shoulder = nsh*fs./newn./newbw;
%          
%            %%% Adjust for upsampling
%            % This will add extra bands at negative frequencies so that the
%            % number of bands is an integer multiple of the frequency upsampling
%            % ratio. This is a bookkeeping measure to ensure that all frequencies
%            % within the desired range are oversampled to the same degree.
%            %
%            %
%            upratio = me.upsampleFx+1;
%            stepsize = newbw/upratio;
%            nstepsize = round(stepsize*newT);
% %            me.offset = me.offset - floor(nsh./(nstepsize+1))*stepsize; % Include any bands that overlap with lower edge so that the lowest band oversampled to the same degree as other bands
% %         foffset = me.offset-newbw*(1+me.shoulder)/2 - me.upsampleFx*stepsize;
%            % This adjusts the high-pass offset so it is an integer multiple of
%            % frequency sampling
%            foffset = floor((me.offset-newbw*(1+me.shoulder)/2 - me.upsampleFx*stepsize)*newT)/newT;
%            noffset = round(foffset*newT );
%            foffset=noffset/newT;
%            me.offset = foffset + newbw*(1+me.shoulder)/2;
%            
%            lowp =  foffset + ceil((me.lowpass -foffset + upratio*stepsize )/(stepsize*upratio))*upratio*stepsize ;% same for highest band
%            
%            
%           
%            me.bands(:,1) = (foffset:stepsize:lowp - stepsize);%-newfs/newn;
%            me.bands(:,2) = me.bands(:,1)+newbw*(1+me.shoulder);%+newfs/newn;
%            nwin =size(me.bands,1);
% 
%            %%% Reshaping matrix. This step includes the initial
%            %%% circular shift.
%            rsmat = noffset +...
%                    repmat((0:nwin-1)*nstepsize,round(winN*(1+me.shoulder)),1) +...
%                    repmat((1:round(winN*(1+me.shoulder)))',1,nwin);% -nsh;
%                          
%            rsmat = mod(rsmat-1,newn)+1;
%            %dcindx = find(rsmat==1);
% 
%              tp = me.taper.make((0:1:nsh-1)/nsh); 
%             invtaper = me.taper.make(1-(0:1:nsh-1)/nsh);
% 
%            switch me.direction
%             %%% Approximate causal or anti-causal filters while
%             %%% preserving summation properties of the tapers.
%             %%% Note that a truly causal filter cannot have a zero frequency response anywhere.   
%                case {'causal','anticausal'}
%                    htptp = hilbert([tp,invtaper]);
%                    htptp = htptp./(abs(htptp)+eps);
%                    if strcmp(me.direction,'causal')
%                       htptp = conj(htptp);  
%                    end
%                    tp = htptp(1:length(tp)).*tp;
%                    invtaper = htptp(length(tp)+(1:length(invtaper))).*invtaper;                       
%                case 'acausal'
%                    % do nothing
%                otherwise
%                    error('Unrecognized filter direction, %s.',me.direction)
%            end
% 
%            Frs = zeros([size(rsmat),ncol]);
%            for  k = 1:ncol
%                 f = F(:,k);
%                Frs(:,:,k) = double(f(rsmat));
%                if nsh>0
%                     Frs(end+(1-nsh:0),1:nwin,k) = diag(sparse(tp))*Frs(end+(1-nsh:0),1:nwin,k);
% 
%                    Frs(1:nsh,1:nwin,k) = diag(sparse(invtaper))*Frs(1:nsh,1:nwin,k); 
%                end
% 
%                 %%% Set all negative frequencies to zero
% %                     Frs(rsmat(:, end-1:end)>ceil((newn+1)/2),end,k)=0;
% %                     Frs(rsmat(:, 1:2)>ceil((newn+1)/2),1,k)=0;
% 
%                 %%% Corrections for DC and Nyquist to preserve power
%                 %%% after zeroing negative frequencies
% %                 Frs(rsmat(:,1)==1,1,k)=Frs(rsmat(:,1)==1,1,k)/sqrt(2); 
%               %  Frs(dcindx + (k-1)*numel(rsmat))=Frs(dcindx + (k-1)*numel(rsmat))/sqrt(2); 
% 
%            end
% 
%            padN = floor((me.fftpad+1)*winN*(1+~me.centerDC)*(1+me.shoulder));
%            me.fftpad = 1/(1+me.shoulder)*padN/winN/(1+~me.centerDC)-1;
%            
%            if padN > round(winN*(1+me.shoulder))
%                Frs(padN,:,:) = 0;
%            elseif padN < round((1+me.shoulder)*winN)
%                padN =round( (1+me.shoulder)*winN);
%                me.fftpad = 0;
%            end
%            
%            if me.centerDC
%                Frs = circshift(Frs,-ceil(winN/2*(1+me.shoulder)));  
%            end
%            
%             me.blrep = ifft(Frs)*sqrt(padN)*sqrt(2)/sqrt(upratio);
%                
%             me.sampling_rate = padN/newT;
%                  
%          
%            me.time = (((1:size(me.blrep,1))-1)*newT./size(me.blrep,1))';
%            w = mean(me.bands,2)';
%             me.blrep(:,w > lowp + me.shoulder*bw,:) = 0;
% %             w(w > me.lowpass + me.shoulder*bw) = [];
%            me.frequency = w;
% 
%            if me.remodphase  
%                %%% This remodulates the signal to obtain the correct phase
%                %%% at a given frequency for the original band. It is necessary 
%                %%% when averaging the complex-valued spectrum over time so 
%                %%% that phase aligns correctly.
%                remodulator = exp(1i*2*pi*me.time*(me.frequency - (~me.centerDC)*me.bandwidth*(1+me.shoulder)/2));
%                me.blrep = me.blrep.*repmat(remodulator,[1 1 size(me.blrep,3)]);
%                
%            end
%            
%         end

        %%%%

%         function [data,fs] = signal(me,columnfilter,doHilbert,return_padded)
%             
%             % db.signal
%             % Reconstruct the signal from its band-limited analytic
%             % representation by applying the inverse DBT. 
%             %
%             % db.signal(bands)
%             %
%             % Reconstruct the signal using only the specified bands. 'bands
%             % may be a numeric or logical index. 
%             % Example: 
%             %        db.signal(db.frequency > 100)
%             % returns the signalusing only bands with center frequency above 100 Hz.
%             %
%             % db.signal(bands,true)
%             %
%             % Returns a complex-valued analytic signal with negative frequencies zero'd. 
%             %
%             
%             
%              upratio = me.upsampleFx+1;
%            
%             if round(me.upsampleFx) ~= me.upsampleFx
%                 error('\n%s does not currently allow signals to be reconstructed from a transform upsampled by a non-integer.',upper(mfilename))
%             end
%             if nargin < 2 || isempty(columnfilter)
%                 mult = diag(sparse(ones(1,upratio)));
%             elseif islogical(columnfilter)
%                 mult = diag(sparse(columnfilter));
%             elseif min(size(columnfilter)) == 1 && min(columnfilter)>=1
%                 
%                mult = diag(sparse( ismember(1:size(me.blrep,2),columnfilter)));
%                
%             elseif min(size(columnfilter)) == 1
%                 mult = diag(sparse(columnfilter));
%             else
%                 mult = columnfilter;
%             end
%             
%             if nargin < 3 || isempty(doHilbert)
%                doHilbert = false; 
%             end
%             if nargin < 4 || isempty(return_padded)
%                return_padded = false; 
%             end
% %             n = me.fullN;
% %             noffset = round((me.offset-me.bandwidth*me.shoulder)./me.fullFS*n);
% %             noffset = round((me.offset - me.bandwidth*(1+me.shoulder)/2)*me.fullN/me.fullFS );
%             ncol = size(me.blrep,3);
%            if me.remodphase  
%                %%% If phase remodulation was applied we need to reverse it.
%                demodulator = exp(-1i*2*pi*me.time*(me.frequency - (~me.centerDC)*me.bandwidth*(1+me.shoulder)/2));
%                me.blrep = me.blrep.*repmat(demodulator,[1 1 size(me.blrep,3)]);
% 
%            end
% 
%             F = zeros(size(me.blrep,1),size(me.blrep,2)/upratio,ncol,upratio);
%             for k = 1:ncol
%                 for upsi = 1:upratio                 
%                     F(:,:,k,upsi) = fft(me.blrep(:,upsi:upratio:end,k) )*mult(upsi:upratio:end,upsi:upratio:end)/sqrt(size(me.blrep,1)/2);         
%                 end            
%             end
%             nsh = round(me.shoulder*me.bandwidth./me.fullFS*me.fullN);
%             padN = length(me.time);
%             winN = round(padN./(1+me.fftpad)/(1+~me.centerDC)/(1+me.shoulder)); 
%             if me.centerDC
%                
%                 F = circshift(F,ceil(winN/2*(1+me.shoulder)));%sqrt(2);
%                 
%             end
%             F = F(1:round(winN*(1+me.shoulder)),:,:,:)*sqrt(2);
%             
%             
%              nnyq = size(F,1);
%             if nsh >0
%                 % Taper is normally defined so that h(k).^2 + h(k+bw).^2 = 1
%                 tp = me.taper.make((0:1:nsh-1)/nsh); 
%                 invtaper = me.taper.make(1-(0:1:nsh-1)/nsh);
%                 switch me.direction
%                 %%% Approximate causal or anti-causal filters while
%                 %%% preserving summation properties of the tapers.
%                    
%                    case {'causal','anticausal'}
%                        htptp = hilbert([tp,invtaper]);
%                        htptp = htptp./(abs(htptp)+eps);
%                        if strcmp(me.direction,'anticausal')
%                           htptp = conj(htptp);  
%                        end
%                        tp = htptp(1:length(tp)).*tp;
%                        invtaper = htptp(length(tp)+(1:length(invtaper))).*invtaper;                       
%                     case 'acausal'
%                         % Do nothing.
%                    otherwise
%                        error('Unrecognized filter direction, %s.',me.direction)
%                 end
%                
%                 for k = 1:ncol
%                     for upsi = 1:upratio
%                         sh = diag(sparse(invtaper))*F(1:nsh,2:end,k,upsi);                
%                         F(nnyq-nsh+1:nnyq,1:end-1,k,upsi) = diag(sparse(tp))*F(nnyq-nsh+1:nnyq,1:end-1,k,upsi)+sh;
%                     end
%                 end
%                  
% %             else
% %                 invtaper = 1;
%             end
%             
% %             if nsh > 0
% %                 %%% Retain the leading edge if there is an offset
% %                F0 = permute(F(1:nsh,1,:),[1 3 2]).*repmat(invtaper',1,ncol); 
% %             else
% %                 F0=zeros(0,ncol);
% %             end
%             
%             F(1:nsh,:,:,:)=[];
%             Ffull = zeros(me.fullN,ncol);
%             switch me.padding
%                 case {'frequency','fft'}
%                     warning('Frequency padding is an obsolete option.')
%             end
%             T = me.fullN./me.fullFS; %Signal duration after padding
%             for k = 1:ncol
%                 
%                for upsi = 1:upratio
%                     f = F(:,:,k,upsi);
% %                Q( 1 + mod(nofs + (1:numel(f))-1,me.fullN),upsi)=f(:);
%                      nofs = round((me.bands(upsi,1) + me.bandwidth*me.shoulder)*T);
%                      %%% Reconstruct signal, averaging over the oversampled
%                      %%% bands.
%                      Ffull( 1 + mod(nofs + (1:numel(f))-1,me.fullN),k) =...
%                         Ffull( 1 + mod(nofs + (1:numel(f))-1,me.fullN),k)...
%                         + f(:) * sqrt(me.fullN)/sqrt(upratio);
%                end
%                if ~mod(me.fullN,2)
%                   Ffull(ceil(me.fullN/2+1),k) = Ffull(ceil(me.fullN/2+1),k)/2;
%                end
%                Ffull(1,k) = real(Ffull(1,k))/2;
%             end
%             Ffull(ceil((me.fullN+1)/2)+1:end,:,:) = 0;
%             
%             Ffull = Ffull./sqrt(2);
%             if doHilbert %%% If second argument is true return analytic signal
%                 data = ifft(full(Ffull));
%             else
%                 data = real(ifft(full(Ffull)));
%             end
%             if ~return_padded
%                 data = data(1:me.Norig,:);
%             end
%             fs = me.fullFS;
%             
%         end
       %%%%
%         function varargout = specgram(me,normalize)
%             if nargin < 2
%                 normalize = 0;
%             end
%             if normalize
%                 fun = @rmbaseline;
%             else
%                 fun = @(x)x.blrep;
%             end
%             S = 20*log10(abs(fun(me)))';
%             t = me.time;
%             w = me.frequency;
%             if nargout == 0
%                imagesc(t,w,S) 
%                axis xy
%             else
%                 varargout(1:3) = {S,t,w};
%             end
%         end
    end
end
    
