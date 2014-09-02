classdef dbt

% The band limited data class for efficient representation of
% band-limited data using the Nyquist sampling theorem. 
% 
% Use:
%
%  B = dbt(X,Fs,BW)
%
%  X  - Signal as a column vector
%  Fs - Sampling frequency
%  BW - Bandwidth of decomposition
% 
%
%  B = dbt(X,Fs,Bw, ['option'], [value])
%  Options:
%
%  	offset   -  offset of the first band from 0hz (default = 0)
%       padding  - 'time': pad signal in the time domain, changing duration (default)
%                  'frequency': pad in the frequency domain, changing sampling rate
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


%     C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

    properties 
        
        blrep = [];  %%% Band-limited analytic representaiton of the signal (time x frequency band)
        bandwidth = [];  %%% bandwidth parameter
        sampling_rate = []; %%% Sampling rate after downsampling
        offset = 0;	%%% Offset if the first band doesn't begin at 0 (a highpass cutoff)
        time = [];     %%% Time points for the rows in blrep
        frequency = [];  %%% Center frequencies for the columns in blrep
        bands =[];     %%% Bands for the columns of blrep 
        fullN = 0;    %%% Length of the signal after padding
        fullFS = 0;   %%% Sampling frequency of the reconstructed signal
        Norig = 0;    %%% Original signal length before padding
        nyqval = 0;   %%% fft value at the nyquist frequency
        shoulder = 1; %%% Degree of frequency overlap between neighboring bands (0 - 1)
        lowpass = []; %%% Lowpass cutoff
%        taperfun =[];
        taper = [];
        padding = 'frequency';
  		 fftpad = 0; % Additional padding to add to the fft as a proportion of unpadded length.
        userdata=[];
        centerDC = false; %%% If false, the pass band is demodulated 
                          %%% so that 0 is the highpass limit, meaning that negative 
                          %%% frequencies contain zero padding.
                          %%% If true, the pass band is centered at 0, so
                          %%% that the signal contains negative
                          %%% frequencies, which halves the number of
                          %%% samples.
                          
%         taper = 'quadratic'; 
    end
    
    
    methods
    
        function me = dbt(varargin)
            
%             padding = 'time';
            i = 4;       
%             me.taperfun = @(x)x;
%            me.taperfun = @(x)(1-cos(x*pi))/2; % Function that defines tapering in the overlapping regions of the window
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
%            me.decim = decim;
           
           bandwindows = 0:bw:fs/2-bw;
           
           fullsig  = varargin{1};
           nyq = bandwindows(end)+bw;
           
           n = size(fullsig,1);
           ncol = size(fullsig,2);
           
           me.Norig = n;
           %Resample signal so that everything factors
           %Keeping signal duration fixed and allowing bandwidth and sampling
           %rate to vary.
           
           T = n./fs; % signal duration
           %%% K and M need to be integers
           
           M = ceil(bw*T);
%            K = ceil(1/2*n./(T*bw) - 1/2*me.offset/bw); 
           K = ceil(1/2*n./M - 1/2*me.offset/bw); 
           newnf = 2*K*M;
           newn = 2*ceil(newnf/2 + me.offset*newnf/(fs-2*me.offset));
           noffset = (newn-newnf)/2;
            
           switch me.padding
               case 'time'
                  T = newn./fs;
                 fullsig(newn) = 0; 
              
                 
               case 'frequency'
                   
                   
           
           end
           
           me.fullN = newn;
           %%% Adjusting bandwidth based on new signal length 
           newbw = 1/2*newnf/(T*K);           
           %%% Likewise for new sampling frequency
           newfs = newn./T;
           me.fullFS = newfs;
           
           if isempty(me.lowpass)
              me.lowpass = newfs/2; 
           end
           
%            nwin = K+2;
           winN = M;
         
%            noffset = floor(me.offset*newn/newfs);
           me.offset = noffset*newfs/newn;
           
        
           F = fft(fullsig./sqrt(n));
           
%            if strcmp(me.padding,'time')
%                F(1,:) = F(1,:) + 1i*F(newn/2+1,:);
%            end
             me.bandwidth = newbw;
           if me.shoulder ~=0
               me.bands(:,1) = [me.offset,me.offset + me.shoulder/2*newbw:newbw:me.lowpass-newbw,me.lowpass- me.shoulder/2*newbw];
               me.bands(2:end-1,2) = me.bands(2:end-1,1)+bw;
               me.bands([1 end],2) = me.bands([1 end],1)+bw*me.shoulder/2;
           else
               me.bands = [me.offset :newbw:me.lowpass-newbw;...
                           (me.offset + newbw) :newbw:me.lowpass ]';
               
           end
           nwin =size(me.bands,1);
           newF = zeros(newn,ncol);           
           oldn = floor(nyq./fs*n)+1;
           newF(1:oldn,:) = F(1:oldn,:);
           
           if me.shoulder == 0
%                nwin = nwin-2;
                nwin = K;
               Frs = reshape(newF(noffset+1:newn/2,:),winN,nwin,ncol);
           else
              
               nsh = min(ceil(me.shoulder*newbw./newfs*newn),winN);
               me.shoulder = nsh*newfs./newn./newbw;
               rsmat = noffset + repmat((0:nwin-3)*(winN),winN+nsh,1) + repmat((1:winN+nsh)',1,nwin-2);
               rsmat = rsmat(:,[1, 1:end,end]);
               tp = me.taper.make((1:1:nsh)/nsh); 
               invtaper = me.taper.make(1-(1:1:nsh)/nsh);
               Frs = zeros([size(rsmat),ncol]);
               for  k = 1:ncol
                    f = newF(:,k);
                   Frs(:,:,k) = double(f(rsmat));
                   Frs(1:nsh,end,k) = diag(sparse(invtaper))*Frs(end-nsh+1:end,end,k);
                   Frs(1:nsh,1,k) = diag(sparse(tp))*Frs(1:nsh,1,k);
                   Frs(nsh+1:end,[1 end] ,k) = 0;
                   %Frs(1:nsh,end ,k) = 0;
                  % Frs(1:nsh,end,k) = diag(sparse(invtaper))*Frs(1:nsh,end,k);
                   
                   %%% Now add the taper
    %                Frs(end+(1-nsh:0),1:nwin-1) = diag(sparse(taper))*Frs(end+(1-nsh:0),1:nwin-1);
                     Frs(end+(1-nsh:0),1:nwin-1,k) = diag(sparse(tp))*Frs(end+(1-nsh:0),1:nwin-1,k);

                   %%% subptract the tapered component from the next band
                   Frs(1:nsh,2:nwin-1,k) = diag(sparse(invtaper))*Frs(1:nsh,2:nwin-1,k); % - Frs(end+(1-nsh:0),1:nwin-1);
                   winN = size(Frs,1);
               end
           end
           padN = floor((me.fftpad+1)*winN*(1+~me.centerDC));
           me.fftpad = padN/winN/(1+~me.centerDC)-1;
           
           if padN > winN
               Frs(padN,:,:) = 0;
           elseif padN < winN
               padN = winN;
               me.fftpad = 0;
           end
%            me.nyqval = newF(newn/2+1,:);
           if me.centerDC
%                Frs = fftshift(Frs,1);
               Frs = circshift(Frs,-ceil(winN/2));
               me.blrep = ifft(Frs)*sqrt(padN);
           else
%                Frs(winN*2,:,:) = 0;
               me.blrep = ifft(Frs)*sqrt(padN);
           end
           
               
            me.sampling_rate = padN/T;
                 
         
%            me.bands(1,1) = me.offset;
%            me.bands(end,2) = me.lowpass;
           me.time = (((1:size(me.blrep,1))-1)*T./size(me.blrep,1))';
%            w = ((0:K-2)+.5 + me.shoulder/2)*newbw + me.offset; 
%            w = [me.shoulder/2*newbw, w, me.lowpass-me.shoulder/2*newbw];
           w = mean(me.bands,2)';
           me.blrep(:,w>me.lowpass,:) = [];
           w(w>me.lowpass) = [];
           me.frequency = w;
           
        end

        %%%%

        function [data,fs] = signal(me,columnfilter,hilbert)
            
            %%% Reconstruct the signal from its band-limited analytic representation
            if nargin < 2 || isempty(columnfilter)
                mult = 1;
            elseif islogical(columnfilter)
                mult = diag(sparse(columnfilter));
            elseif min(size(columnfilter)) == 1
                
               mult = diag(sparse( ismember(1:size(me.blrep,2),columnfilter)));
               
            else
                mult = columnfilter;
            end
            
            if nargin < 3 || isempty(hilbert)
               hilbert = false; 
            end
            
            n = me.fullN;
            noffset = round(me.offset./me.fullFS*n);
            ncol = size(me.blrep,3);
            F = zeros(size(me.blrep));
            for k = 1:ncol
                F(:,:,k) = fft(me.blrep(:,:,k) )*mult/sqrt(size(me.blrep,1)/2);         
            end
            nsh = round(me.shoulder*me.bandwidth./me.fullFS*me.fullN);
            padN = length(me.time);
            winN = round(padN./(1+me.fftpad)/(1+~me.centerDC)); 
            if me.centerDC
%                 nnyq = size(F,1);
                F = circshift(F,ceil(winN/2));%sqrt(2);
                
%                 F = ifftshift(F,1)*sqrt(2);
%             else
               % nnyq = ceil((size(F,1)+1)/2);
            end
            F = F(1:winN,:,:)*sqrt(2);

             nnyq = size(F,1);
            if nsh >0
                % Taper is normally defined so that h(k).^2 + h(k+bw).^2 = 1
                tp = me.taper.make((1:1:nsh)/nsh); 
                invtaper = me.taper.make(1-(1:1:nsh)/nsh);
                for k = 1:ncol
                    sh = diag(sparse(tp))*F(nnyq-nsh+1:nnyq,1:end-1,k);
                    sh(:,1) = diag(sparse(tp))*F(1:nsh,1,k);
                    %sh(:,end) = (diag(sparse(invtaper))*F(nnyq-nsh+1:nnyq,end,k)+sh(:,end));
                    
                    F(1:nsh,2:end,k) = diag(sparse(invtaper))*F(1:nsh,2:end,k)+sh;
                end
                    F(:,1,:) = [];
                 
            end
            F(nnyq-nsh+1:end,:,:) = [];
            Ffull = zeros(me.fullN,ncol);
            switch me.padding
                case {'frequency','fft'}
                    for k = 1:ncol
                      f = F(:,:,k);

                       Ffull(noffset+(1:numel(f)),k) = f(:)*sqrt(me.Norig);
                       Ffull(floor(me.Norig/2+.5)+1,k) = Ffull(floor(me.Norig/2+.5)+1,k)/2; 
                       Ffull(ceil(me.Norig/2+.5)+1:me.Norig,k) = 0; 

                        Ffull(1,k) = real(Ffull(1,k))/2;
                    end
                    
                    Ffull(me.Norig+1:end,:) = []; 

%                     Ffull(ceil(me.Norig/2)+1) = imag(Ffull(1))/2;
%                     Ffull(1) = real(Ffull(1))/2;
                 case 'time'
                    for k = 1:ncol
                        f = F(:,:,k);
                        Ffull(noffset+(1:numel(f)),k) = f(:)*sqrt(me.fullN);

                        Ffull(me.fullN/2+1,k) = imag(Ffull(1,k))/2;
                        Ffull(1,k) = real(Ffull(1,k))/2;
%                         Ffull(me.fullN,k) = 0;
                    end
            end

            if hilbert
                data = ifft(full(Ffull));
            else
                data = real(ifft(full(Ffull)));
            end
            data = data(1:me.Norig,:);
            fs = me.fullFS;
            
        end
       
        function varargout = specgram(me,normalize)
            if nargin < 2
                normalize = 1;
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
            else
                varargout(1:3) = {S,t,w};
            end
        end
    end
end
    
