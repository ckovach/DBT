function [coh,csp,w,tt,dbs,trf,Pperm] =dbtcoh(x,y,varargin)


% coh =dbtcoh(x,y,fs,bw)
% Efficiently compute coherences values with dbt transforms. x and y are
% input signals arranged time x channel. fs is sampling rate and bw is the
% bandwidth within which to compute coherence.
%
% If x and y are dbt objects (see DBT), then fs and bw can be ommitted.
%
% [coh,csp,w] =dbtcoh(x,y,fs,bw)
% Returns the cross-spectral matrix and vector of center frequencies for each band, w.
%
% [coh,csp,w,tt] =dbtcoh(x,y,fs,bw,'trigger',trig,'timerange',trange)
%
% Computes event-related coherence by averaging over windows specified by
% trigger and timerange. Trigger is a vector of event times and time range
% is a vector, trange = [startt endt], specificying window onset and end
% relative to events in trig. tt is a vector sample times for the window.
%
% See also DBT

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

i = 1;
keep_time=[];
trigger = [];
timerange=[];
tt=0;
subtract_mean = true; % Subtract the average as with an ordinary correlation.
                       % This generally should make little difference on
                       % simple coherence analyses but may be important for
                       % event-related PAC.
permtest = nargout>6; % Run a permutation test if true
nperm = 1000; % number of permutations

dbtargs = {'remodphase',subtract_mean};


if nargin>3 && isnumeric(varargin{2}) && isnumeric(y)
    fs = varargin{1};
    bw = varargin{2};
    varargin(1:2) =[];
elseif nargin > 2 && isscalar(y) && isnumeric(y)
    fs = y;
    bw = varargin{1};
    varargin(1) = [];
end
                       
while i <= length(varargin)
  
    
   switch varargin{i}
       
       case 'keep time'
           keep_time = varargin{i+1};
           i = i+1;
       case 'trigger'
           trigger = varargin{i+1};
           i = i+1;
       case 'timerange'
           timerange = varargin{i+1};
           i = i+1;
       case {'subtract mean','center'}
           subtract_mean = varargin{i+1};
           dbtargs = [dbtargs,{'remodphase',subtract_mean}]; %#ok<*AGROW>
           i = i+1;
       otherwise
           dbtargs = [dbtargs,varargin(i:i+1)]; 
           i=i+1;
   end
   
   i=i+1;
end


if ~isa(x,'dbt')
 
    nx = size(x,2);   
    
else
    dbx = x;
    nx = size(dbx.blrep,3);
end

 dby = [];
 
if ~isempty(y) && (~isscalar(y) || isa(y,'dbt')) 
    if ~isa(y,'dbt')
        dby = dbt(y,fs,bw,dbtargs{:});
        ny = size(y,2);
    else
        dby = y;
        ny = size(dby.blrep,3);
    end
else
    ny = nx;    
end

if ~isa(x,'dbt')
    dbx = dbt(x,fs,bw,dbtargs{:});
end

if isempty(keep_time)
    keepT = true(size(dbx.time));
else
    keepT = resampi(keep_time,fs,dbx.sampling_rate,'linear')>.5;
end

w = dbx.frequency;
csp = zeros(nx,ny,length(w));
coh = csp;
if nargout > 4
    trf = csp;
end

if ~isempty(timerange)
   [AX,tt] =  choptf(timerange,trigger,dbx);
   if ~isempty(y)
       [AY,tt] =  choptf(timerange,trigger,dby);
   end
else
    AX = dbx.blrep(keepT,:,:);
    AY = dby.blrep(keepT,:,:);
end


for i = 1:length(dbx.frequency)
    
    
    for t = 1:length(tt)
        if isempty(timerange)
            blx = squeeze(AX(:,i,:));
            if isempty(y)
                bly = blx;
            else
                bly = squeeze(AY(:,i,:));
            end
        else
            blx = squeeze(AX(t,i,:,:));
            if isempty(y)
                bly = blx;
            else
                bly = squeeze(AY(t,i,:,:));
            end

        end
        csp(:,:,i,t) = blx'*bly;
        
        if permtest
            nt = size(blx,1);
            permN = 0;
           for k = 1:nperm
               rp = randperm(nt);
               permN = permN + (abs(blx(rp,:)'*bly)>abs(csp(:,:,i,t))); %%% Permutes over time intervals.
           end
            Pperm(:,:,i,t) = permN./nperm;
        else
            Pperm = nan;
        end
        
        if subtract_mean
                 csp(:,:,i,t) = csp(:,:,i,t) - sum(blx)'*mean(bly);
        end
            
        if isempty(y)
             coh(:,:,i,t) = diag(diag(csp(:,:,i,t).^-.5))*csp(:,:,i,t)*diag(diag(csp(:,:,i,t)).^-.5);
        else          
           coh(:,:,i,t) = diag(sum(abs(blx).^2).^-.5)*csp(:,:,i,t)*diag(sum(abs(bly).^2).^-.5);
        end

        if nargout > 4
%             cblx = blx-repmat(mean(blx),size(blx,1),1);
%             cbly = bly-repmat(mean(bly),size(bly,1),1);
%            trf(:,:,i,t) = diag(sum(abs(blx).^2))\(blx'*bly);
           trf(:,:,i,t) = diag(sum(abs(bly).^2))\(blx'*bly);
        end

    end
end

    
if nargout > 4
    if isequal(x,y)
        dbs = dbx;
    else
        dbs = [dbx,dby];
    end
end



