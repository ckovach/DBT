function [xypac,amp,rmph] = dbtpac(X,Y,fs,varargin)

% function dbtpac(x,y,fs,varargin)
% Efficiently computes phase-amplitude coupling with the dbt transform
% Phase is computed from x and amplitude from y. Because the amplitdue 
% component is necessarily sampled at a higher bandwidth, the result is
% averaged over intervening steps. 
%
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

phasebw = 4;
ampbw =40; 
phaserange = [0.01 24];
amprange = [40 300];
phasefs=fs;
coherence = false; % compute coherence instead of PAC if true
smwin = 0;
anglehist = false;
partial =  false;
exclude_time = false;
i = 1;
while i < length(varargin)
    
   switch varargin{i}
       
       case 'phasebw'
           phasebw = varargin{i+1};
           i = i+1;
       
       case 'ampbw'
           ampbw = varargin{i+1};
           i = i+1;
       case 'coherence'
           coherence = varargin{i+1};
           i = i+1;
       case 'phase range'
           phaserange = varargin{i+1};
           i = i+1;
       case 'amp range'
           amprange = varargin{i+1};
           i = i+1;
       case 'phasefs'
           phasefs= varargin{i+1};
           i = i+1;
       case 'smwin'
           smwin= varargin{i+1};
           i = i+1;
         case 'partial'
           partial= varargin{i+1};
           i = i+1;
         case 'exclude time'
           exclude_time= varargin{i+1};
           i = i+1;
       otherwise
           error('Unrecognized keyword %i',varargin{i+1})
   end
   i = i+1;
end


    
ampargs = {fs,ampbw,'padding','frequency','lowpass',amprange(2),'offset',amprange(1)};
fprintf('\nAmplitude: %4i',0)
% fprintf('\nAmplitude %4i',0)
for i = 1:size(Y,2)
    fprintf('\b\b\b\b%4i',i)
    y = Y(:,i);
    dbamp = dbt(y,fs,ampbw,'padding','frequency','lowpass',amprange(2),'offset',amprange(1));
%     dbampdb = dbt(abs(dbamp.blrep),dbamp.sampling_rate,phasebw,'padding','frequency','lowpass',phaserange(2),'offset',phaserange(1));
    a = zscore(bpfilt(abs(dbamp.blrep),[dbamp.sampling_rate phaserange ]));  
    w = (0:size(a,1)-1)/size(a,1)*dbamp.sampling_rate;
    Fa = fft(a);
    Fa(w>phaserange(2)*2,:) = [];
%     amp(:,i,:) = real(ifft(Fa));
    amp(:,i,:) = ifft(Fa);
    if i==1 && size(Y,2)>1
        amp(end,size(Y,2),end) = 0;
    end
        
end
ampfs = size(Fa,1)./size(a,1)*dbamp.sampling_rate;
ampt = (0:size(amp,1)-1)./ampfs;

%         fprintf('\b\b\b\b\b\b\b\b\b\b\b')
%         fprintf('%0.4i<->%0.4i',i,k)

    fprintf('\nPhase: %4i',0)
phargs = {phasefs,phasebw,'padding','frequency','lowpass',phaserange(2),'offset',phaserange(1)};
for k = 1:size(X,2);        
       fprintf('\b\b\b\b%4i',k)

        x = X(:,k);
        dbph = dbt(x,phargs{:});


        %Interpolate to amplitude sampling rate with sinc interpolation (i.e. padding the fft).
        ph(:,k,:) = dbph.blrep;
        if k==1 && size(X,2)>1
            ph(end,size(X,2),end) = 0;
        end 

%          camp = 0;
end
if islogical(smwin)&&smwin 
   smwin =  1/dbph.sampling_rate;
end

nx = size(X,2);
ny= size(Y,2);
nphbands = length(dbph.frequency);
nsamp = length(ampt);
nampbands = length(dbamp.frequency);


%%% Match phase blrep to amplitude blrep sampling rate with sinc interp
%%% (i.e. pad or trim the fft)
F = fft(ph);
F(nsamp,:,:,:) = 0;
iph = ifft(F);

if ~coherence
   iph = iph./abs(iph);
end

%%% Construct remodulator to adjust phase to a correct value 
%%% at the sampled time
rmf = dbph.bands(:,1)-dbph.shoulder/2*dbph.bandwidth;
remodulator = repmat(permute(exp(2*pi*1i*rmf*ampt)',[1 3 2]), [1 nx  1  ]);

if ~isscalar(exclude_time)
   ext = (0:length(exclude_time))/fs;
   gett = interp1(ext,double(exclude_time),ampt,'nearest')>0;
else
    gett = true(size(ampt));
end
%%% remodulate phase and reshape matrices
rmph = iph.*conj(remodulator);
% amp = zscore(amp);
for i = 1:nampbands
    fprintf('%4i',0)
    
    for k = 1:nphbands
%         fprintf('%4i',k)
        if partial
            a = amp(:,:,i);
%             p = rmph(:,:,k);
            b =sum(a.*conj(rmph(:,:,k)))./sum(abs(rmph(:,:,k)).^2);
            
            a = a-rmph(:,:,k).*repmat(b,size(a,1),1);
        else
            a = amp(:,:,i);
        end
            
        PAC(:,:,i,k) = corr(a(gett,:,:),rmph(gett,:,k));

        if i==1 && ny > 1
            PAC(:,:,nampbands,nphbands)= 0;
        end
%         fprintf('\b\b\b\b')
    end
    
    fprintf('\b\b\b\b')
    
end
% CM = rmph.*real(amp);
% PAC(:,:,i,:) = squeeze(sum(CM)./sqrt(sum(abs(amp).^2).*sum(abs(rmph).^2)));

if smwin > 0
   %%% If cmwin is greater than zero, include a time-varying comodulogram in the output 
   %%% To construct this, note that significant PAC implies that the
   %%% phase component of one signal demodulates the other, hence we can
   %%% low-pass filter to obtain a time-varying version of this. Here we
   %%% will simply downsample to match dbph.
        smwinN = round(smwin*dbamp.sampling_rate);
           smw = zeros(size(amp,1),1);
           smw(1:smwinN,:) = hanning(smwinN);
           smw = circshift(smw,-floor(smwinN/2));
           smw = smw./sum(smw);

           fsmw = fft(smw);
   for k = 1:size(amp,3)
       for  kk = 1:size(rmph,3)
           CM = repmat(amp(:,:,k),[1 1 size(rmph,2)]).*repmat(permute(conj(rmph(:,:,kk)),[1 3 2]),[1 size(amp,2) 1]);
 
           
   
           rp = size(CM);rp(1)=1;
           FCM= fft(CM).*conj(repmat(fsmw,rp));
           w = (0:nsamp-1)/nsamp*dbamp.sampling_rate;
           FCM = FCM(1:length(dbph.time),:,:,:,:,:);

           Famp = fft(abs(amp(:,:,k)).^2).*conj(repmat(fsmw,1,size(amp,2)));
           Famp = Famp(1:length(dbph.time),:,:,:,:,:);

           CMrs(:,:,:,k,kk) = ifft(FCM)./repmat(sqrt(abs(ifft(Famp))),[1 1 size(rmph,2)]);
           if k==1 && kk == 1 && size(amp,3)>1 && size(rmph,3)>1
               CMrs(:,:,:,size(amp,3),size(rmph,3)) = 0;
           end
       end
   end
   cmtime = ampt;
else
    CMrs = [];
    cmtime = [];
end
fprintf('\b\b\b\b')

xypac.PAC = PAC;
% xypac.amp = amp;
% xypac.ph = rmph;
xypac.PACgram =CMrs;
xypac.ampfreq = dbamp.frequency';
xypac.phfreq = dbph.frequency;
xypac.fs = fs;
xypac.cmtime= cmtime;
xypac.args.in = varargin;
xypac.args.dbamp = ampargs;
xypac.args.dbphase = phargs;


