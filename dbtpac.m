function [xypac,dbamp,dbph] = dbtpac(X,Y,fs,varargin)

% out = dbtpac(x,y,fs,'phasebw',phbw,'ampbw',abw,varargin)
% 
% Efficiently computes phase-amplitude coupling, measured as LFP-amplitude 
% coherence values, using the dbt transform.
%
% First, band-limited amplitude envelopes are computed from Y with the DBT at
% bandwidth abw (default 40 Hz), then X is downsampled to the amplitude-band
% sampling rate and amplitude-LFP coherence is calculated with DBTCOH at
% bandwidth phasebw (defaults to 50 / T, where T is the duraion 
% of the signal). 
% 
% The output, out, is a struct with phase-amplitude coherence in the field out.PAC. 
% out.PAC is a 4-d matrix with the result arranged as channel(phase) x
% channel(amplitude) x amplitude bands x phase bands.
% 
% See also DBTCOH, DBT

% C. Kovach 2014

phasebw = [];
ampbw =40; 
phaserange = [];
amprange = [0 300];

% coherence = false; % compute coherence instead of PAC if true
% smwin = 0;
% anglehist = false;
% partial =  false;
% exclude_time = false;
% timerange = [];
% trigger= [];
keep_time = [];
phargs ={};
get_trf = false;
get_csd = false;
trf = [];
csd =[];
cohargs = {};
dbtargs = {};
do_permtest = false;
% if nargin <3 && isa(X,'dbt')
%     fs = X.fullFS;
% end
% phasefs=fs;

i = 1;
while i < length(varargin)
    
   switch varargin{i}
       
       case 'phasebw'
           phasebw = varargin{i+1};
           i = i+1;
       
       case 'ampbw'
           ampbw = varargin{i+1};
           i = i+1;
%        case 'coherence'
%            coherence = varargin{i+1};
%            i = i+1;
       case 'phase range'
           phaserange = varargin{i+1};
           i = i+1;
       case 'amp range'
           amprange = varargin{i+1};
           i = i+1;
%        case 'phasefs'
%            phasefs= varargin{i+1};
%            i = i+1;
%        case 'smwin'
%            smwin= varargin{i+1};
%            i = i+1;
%          case 'partial'
%            partial= varargin{i+1};
%            i = i+1;
         case 'keep time'
           keep_time= varargin{i+1};
           i = i+1;
         case {'transfer function','trf'}
           get_trf= varargin{i+1};
           i = i+1;
         case {'csd','cross spectrum'}
           get_csd= varargin{i+1};
           i = i+1;
        case {'trigger','timerange'}
           phargs = [phargs,varargin(i:i+1)];
           i = i+1;
        case {'center','subtract mean'}
           cohargs = [cohargs,varargin(i:i+1)];
           i = i+1;
        case {'dbtargs'}
           dbtargs = [dbtargs,varargin{i+1}]; %#ok<*AGROW>
           i = i+1;
       case 'do_permtest'
           do_permtest = varargin{i+1};
           i = i+1;
       otherwise
           error('Unrecognized keyword %s',varargin{i})
   end
   i = i+1;
end


    
ampargs = {fs,ampbw,'padding','time','lowpass',min(amprange(2),fs/2),'offset',amprange(1),dbtargs{:}};
dbamp = dbt(Y,ampargs{:}); %%% DBT from which band-limited amplitude will be obtained. 
ampfs = dbamp.sampling_rate;

if isempty(phaserange)
   phaserange = [0 dbamp.sampling_rate/2]; 
end

if isempty(phasebw)
    %%% If no other bandwidth is specified default to a reasonable value
    %%% based on the amplitude bandwidth.
   phasebw = ampbw/length(dbamp.time)*100; 
end

if isempty(keep_time)
    keepT = dbamp.time<=dbamp.Norig./dbamp.fullFS;
else
    keepT = resampi(keep_time,fs,ampfs,'linear')>.5;
end
% Downsample data to match the DBT sampling rate for amplitude. 
Xrs = resampi(X,fs,ampfs,'fft');
Xrs(end+1:length(dbamp.time),:,:) =0;
Pperm =[];
phargs = {phasebw,'padding','time','lowpass',min(phaserange(2),fs/2),'offset',phaserange(1),'keep time',keepT,phargs{:},dbtargs{:}}; %#ok<*CCAT>
fprintf('\nBand: %4i',0)
for k = 1:length(dbamp.frequency)   
    fprintf('\b\b\b\b%4i',k)
    if get_trf
        [PAC(:,:,k,:,:),c,phfreq,tt,dbph,trf(:,:,k,:,:)] = dbtcoh(Xrs,squeeze(abs(dbamp.blrep(:,k,:))),ampfs,phargs{:},cohargs{:});
    elseif do_permtest
         [PAC(:,:,k,:,:),c,phfreq,tt,dbph,trf(:,:,k,:,:),Pperm(:,:,k,:,:)] = dbtcoh(Xrs,squeeze(abs(dbamp.blrep(:,k,:))),ampfs,phargs{:},cohargs{:});        
    else
        [PAC(:,:,k,:,:),c,phfreq,tt,dbph] = dbtcoh(Xrs,squeeze(abs(dbamp.blrep(:,k,:))),ampfs,phargs{:},cohargs{:});        
    end
    if get_csd
        csd(:,:,k,:,:) = c;
    end
end

fprintf('\b\b\b\b')

xypac.PAC = PAC;
% xypac.amp = amp;
% xypac.ph = rmph;
% xypac.PACgram =CMrs;
xypac.ampfreq = dbamp.frequency';
xypac.phfreq = phfreq;
xypac.fs = fs;
% xypac.cmtime= cmtime;
xypac.args.in = varargin;
xypac.args.dbamp = ampargs;
xypac.args.dbphase = phargs;
xypac.tt = tt;
xypac.trf = trf;
xypac.csd = csd;
xypac.Pperm= Pperm;
