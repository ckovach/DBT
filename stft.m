classdef stft

% STFT signal representation. 
% 
% Use:
%
%  B = stft(X,Fs,TW)
%
%  X  - Signal as a column vector
%  Fs - Sampling frequency
%  TW - Time window width 
% 
%
%  B = stft(X,Fs,TW, ['option'], [value])
%  Options:
%
%  	    offset   -  offset of the first band from 0hz (default = 0)
%       padding  - 'time': pad signal in the time domain, changing duration (default)
%                  'frequency': pad in the frequency domain, changing sampling rate
%       shoulder - (0 - 1) degree of overlap between neighboring bands (default = 0)           
%
%  The overlapping portions of the bands are windowed with a linear taper:
%	     ____   __________   ________   _________
%	         \ /	      \ /        \ /   
% .  .  .	  X            X          X           . . . 
% 	         / \	      / \        / \
%            |-----------|           |-|
%    	               TW              shoulder
%      


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
        timewindow = [];  %%% bandwidth parameter
        sampling_rate = []; %%% Sampling rate after downsampling
        time = [];     %%% Time points for the rows in blrep
        frequency = [];  %%% Center frequencies for the columns in blrep
        windows =[];     %%% Time ranges for the columns of blrep 
        fullN = 0;    %%% Length of the signal after padding
        fullFS = 0;   %%% Sampling frequency of the reconstructed signal
        Norig = 0;    %%% Original signal length before padding
        nyqval = 0;   %%% fft value at the nyquist frequency
        shoulder = 0; %%% Degree of frequency overlap between neighboring bands (0 - 1)
        lowpass = []; %%% Lowpass cutoff
%        taperfun =[];
        taper = [];
          twtol = 1e-8;    % 
        upsampleTx = 0;
        fftpad = 0;
        times =[];
%         taper = 'quadratic'; 
    end
    
    
    methods
    
        function me = stft(varargin)
            
            padding = 'time';
            i = 4;       
%             me.taperfun = @(x)x;
%            me.taperfun = @(x)(1-cos(x*pi))/2; % Function that defines tapering in the overlapping regions of the window
            me.taper = taper;
            
           while i < length(varargin)
              switch lower(varargin{i})
                  
                  case 'band'
                      band = varargin{i+1};
                      me.offset = band(1);
                      me.lowpass = band(end);                      
                      i = i+1;
                  case 'padding'
                      padding = varargin{i+1};
                      i = i+1;
                  case 'shoulder'
                      me.shoulder = varargin{i+1};
                      i = i+1;
                  case 'lowpass'
                      me.lowpass = varargin{i+1};
                      i = i+1;
                 case 'taper'  %Define a taper function. Default is (1-cos(x))/2
                      me.taperfun = varargin{i+1};
                      i = i+1;
                case {'fftpad','upsamplefx'} % Pad by proportion of window size
                      me.fftpad = varargin{i+1};
                      i = i+1;
                 case 'upsampletx' %Upsample time by proportion + 1 
                      
                      nx =  varargin{i+1};
                    
                      if nx ~= round(nx)
                        warning('%s will not currently allow signals to be reconstructed with non-integer upsampling. Integer upsampling is recommended.',upper(mfilename))
                      end
            
                      me.upsampleFx = nx;
                      i = i+1;                      
                                            
                  
                 case 'upsample' %Upsample time by proportion with fft padding and frequency by increasing overlap
                      nx =  varargin{i+1}-1;
                      me.fftpad = nx;
                      me.upsampleTx = nx;
                      if nx ~= round(nx)
                          warning('%s will not currently allow signals to be reconstructed with non-integer upsampling. Integer upsampling is recommended.',upper(mfilename))
                      end
                      i = i+1;
                  otherwise
                     error('Unrecognized keyword %s',varargin{i})
              end
              i = i+1;
           end
               
           if isempty(varargin)
               return
           end
           
        
           fullsig  = varargin{1};
           fs  = varargin{2};
           tw =  varargin{3};

           n = size(fullsig,1);
           T = n./fs; % signal duration
%            me.decim = decim;
           
            stepsize = tw/(me.upsampleTx +1);           
        %    [nm,den] = rat(T./stepsize,me.twtol);
     
           nchan = size(fullsig,2);
           
           me.Norig = n;
           
%            timewindows = 0:tw-me.shoulder:T-tw;
        
           
           nstepsize = ceil(stepsize*fs);
           stepsize = nstepsize./fs;
           newT = ceil(T./stepsize)*stepsize;
           newn = newT.*fs;
           newtw = stepsize.*(me.upsampleTx+1);
                    
           winN = round(fs*newtw);
           
           
           nsh = round(me.shoulder*tw*fs);
         
           upratio = me.upsampleTx+1;
           stepsize = newtw/upratio;
           nstepsize = round(stepsize.*fs);
                 
      
           
           me.fullN = newn;
           %%% Likewise for new sampling frequency
           newfs = newn./newT;
           me.fullFS = newfs;
           
           
           me.times(:,1) = (0:stepsize:newT - stepsize);%-newfs/newn;
           me.times(:,2) = me.times(:,1)+newtw*(1+me.shoulder);%+newfs/newn;
           nwin =size(me.times,1);

           %%% Reshaping matrix. This step includes the initial
           %%% circular shift.
           rsmat =  repmat((0:nwin-1)*nstepsize,round(winN*(1+me.shoulder)),1) +...
                   repmat((1:round(winN*(1+me.shoulder)))',1,nwin);% -nsh;
           
           fullsig(end+1:max(rsmat(:)),:,:)=0;           
           %dcindx = find(rsmat==1);

             tp = me.taper.make((0:1:nsh-1)/nsh); 
            invtaper = me.taper.make(1-(0:1:nsh-1)/nsh);
        
           
           for  k = 1:nchan
                x = fullsig(:,k);
               Xrs(:,:,k) = double(x(rsmat));
               if nsh>0
                    Xrs(end+(1-nsh:0),1:nwin,k) = diag(sparse(tp))*Xrs(end+(1-nsh:0),1:nwin,k);

                   Xrs(1:nsh,1:nwin,k) = diag(sparse(invtaper))*Xrs(1:nsh,1:nwin,k); 
               end
              
           end
           padN = floor((me.fftpad+1)*winN*(1+me.shoulder));
           me.fftpad = 1/(1+me.shoulder)*padN/winN-1;
      
           Xrs(end+1:padN,:,:) =0; 
           Frs = fft(Xrs);
           me.blrep = permute(2*Frs(1:end/2,:,:),[2 1 3]);
           me.blrep(:,1,:) = Frs(1,:,:);
           
           me.sampling_rate = nwin/newT;
           
                              
           me.timewindow = newtw;
           me.windows = [0:newtw:me.lowpass-newtw;(newtw:newtw:newT)]';
           me.time = mean(me.times,2);
           w = (0:padN-1)./padN*newfs/2;          
           me.frequency = w;
           
        end

        %%%%

        function [data,fs] = signal(me,rowfilter,hilbert)
            
            %%% Reconstruct the signal from its band-limited analytic representation
            if nargin < 2 || isempty(rowfilter)
                mult = 1;
            elseif islogical(rowfilter)
                mult = diag(sparse(rowfilter));
            elseif min(size(rowfilter)) == 1
                
               mult = diag(sparse( ismember(1:size(me.blrep',1),rowfilter)));
               
            else
                mult = rowfilter;
            end
            
            if nargin < 3 || isempty(hilbert)
               hilbert = false; 
            end
            
            F = mult*me.blrep';
            F(2*end,:)=0;
            n = me.fullN;
            Xrs = ifft(F); 
            Xrs = Xrs(1:end/2,:);
            nsh = round(me.shoulder*me.timewindow.*me.fullFS);
            
            % Taper is normally defined so that h(k).^2 + h(k+bw).^2 = 1
            tp = me.taper.make((1:1:nsh)/nsh); 
            invtaper = me.taper.make(1-(1:1:nsh)/nsh);

		    sh = diag(sparse(tp))*Xrs(end-nsh+1:end,1:end-1);
            Xrs(1:nsh,2:end) = diag(sparse(invtaper))*Xrs(1:nsh,2:end)+sh;
            Xrs(end-nsh+1:end,:) =[];
            
                data = Xrs(1:me.Norig);
            if ~hilbert
                data = real(data);
            end
            fs = me.fullFS;
            
        end
        function delete(me)
           me.blrep = []; 
        end
       
    end
end
    
