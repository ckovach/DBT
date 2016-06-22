function [coh,csp,w,tt,dbs,trf,bias,Pperm] =dbtcoh(x,y,varargin)


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
event_type = [];
tt=0;
rm0lag = false; %remove zero phase lag component through a frequency domain highpass filter
                % If rm0lag = true, then subtracts the mean from the cross
                % spectrum over frequencies.
                % If rm0lag is a real scalar, then rm0lag applies a
                % 'highpass' frequency domain filter, removing the
                % contribution within a time lag of the given value.
subtract_mean = true; % Subtract the average as with an ordinary correlation.
                       % This generally should make little difference on
                       % simple coherence analyses but may be important for
                       % event-related PAC.
% permtest = nargout>7; % Run a permutation test if true
nperm = 0; % number of permutations

dbtargs = {'remodphase',subtract_mean};

gpuEnable = false;

% aPLV = false; %Compute amplitude-weighted phase locking rather than coherence.
type  = 'coh';
resolve_time = true;

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
       case 'type'
           type = varargin{i+1};
           i = i+1;
       case 'permtest'
           nperm = varargin{i+1};
           i = i+1;
       case 'conditions'
           event_type = varargin{i+1};
           i = i+1;
       case 'resolve time'
           resolve_time = varargin{i+1};
           i = i+1;
       case 'gpu'
           gpuEnable = varargin{i+1};
           i = i+1;
       otherwise
           dbtargs = [dbtargs,varargin(i:i+1)]; 
           i=i+1;
   end
   
   i=i+1;
end

gpuEnable = gpuEnable && gpuDeviceCount >0;

switch lower(type) 
    case {'blplv','blpl'}
        dbtargs = [dbtargs,{'centerDC',false,'remodphase',false}];
    case {'icplv'}
        dbtargs = [dbtargs,{'remodphase',true}];
        
end

if ~isa(x,'dbt')
 
    nx = size(x,2);   
    
else
    dbx = x;
    nx = size(dbx.blrep,3);
    gpuEnable = gpuEnable || dbx.gpuEnable; % Enable if dbt object is already GPU enabled. 
    
end

if gpuEnable
   dbtargs = [dbtargs,{'gpu',true}]; 
   gpuarg = {'gpuArray'};
else
    gpuarg = {};
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

switch lower(type) 
    case {'blplv','blpl'}
        dbx = blphase(dbx);
        dby = blphase(dby);
end

if isempty(keep_time)
   % keepT = true(size(dbx.time));
   keepT = dbx.time < dbx.Norig/dbx.fullFS;

elseif isa(x,'dbt')
    keepT = keep_time;
else
    keepT = resampi(keep_time,fs,dbx.sampling_rate,'linear')>.5;
end

w = dbx.frequency;
csp = zeros(nx,ny,length(w),gpuarg{:});

if nargout > 4
    trf = csp;
end

if ~isempty(timerange) 
   [AX,tt] =  choptf(timerange,trigger,dbx);
   if ~isempty(y)
       [AY,tt] =  choptf(timerange,trigger,dby);
   end
else
    unqev =0;
    AX = dbx.blrep(keepT,:,:);
    AY = dby.blrep(keepT,:,:);
end

if isempty(event_type)
    event_type = zeros(size(tt));
end
[unqev,~,unqi] = unique(event_type);
 
% bias= zeros(nx,ny,length(w),gpuarg{:});

coh = zeros([nx,ny,length(w),length(tt),length(unqev),nperm+1],gpuarg{:});
if nargout > 6
    bias = zeros([nx,ny,length(w),length(tt),length(unqev),nperm+1],gpuarg{:});
end
for permi = 1:nperm+1
     if  permi>1
             nt = size(AX,1);
            rp = randperm(nt);  
            AX = AX(rp,:,:);
      end
for i = 1:length(dbx.frequency)
    
    for t = 1:length(tt)
       for k = 1:length(unqev)
        if isempty(timerange)
            blx = squeeze(AX(:,i,:));
            if isempty(y)
                bly = blx;
            else
                bly = squeeze(AY(:,i,:));
            end
        else
            if resolve_time
            
             blx = squeeze(AX(t,i,event_type == unqev(k),:));

                 if isempty(y)
                    bly = blx;
                else
    %                 bly = squeeze(AY(:,i,:));
                  bly = squeeze(AY(t,i,event_type == unqev(k),:));
                 end
            else
                  blx = reshape(AX(:,i,event_type == unqev(k),:),[size(AX,1)*size(AX,3), size(AX,4) 1]) ;
                 if isempty(y)
                    bly = blx;
                 else
                    bly = reshape(AY(:,i,event_type == unqev(k),:),[size(AY,1)*size(AY,3), size(AY,4) 1]) ;

                 end
            end
                       
        end
        
        switch lower(type)
            case {'blplv','blpl','plv'}
                blx = blx./abs(blx);
                    bly = bly./abs(bly);
            case {'icplv'}
                blx = itercent(blx./abs(blx));
                bly = itercent(bly./abs(bly));
                
        end
        csp(:,:,i,t,k) = blx'*bly;
        
       
        
        if subtract_mean
                 csp(:,:,i,t,k) = csp(:,:,i,t,k) - sum(blx)'*mean(bly);
        end
            
        if isempty(y)
             coh(:,:,i,t,k,permi) = diag(diag(csp(:,:,i,t,k).^-.5))*csp(:,:,i,t,k)*diag(diag(csp(:,:,i,t,k)).^-.5);
        else     
            switch lower(type)
                case {'awplv','aplv'}
                     coh(:,:,i,t,k,permi) = csp(:,:,i,t,k)./(abs(blx)'*abs(bly));
                    % bias(:,:,i,t,k) = sqrt(sum((abs(blx).*abs(bly)).^2)./sum(abs(blx)'.*abs(bly)).^2);
                   if nargout > 6
                       bias(:,:,i,t,k,permi) = sqrt((abs(blx))'.^2*(abs(bly)).^2)./(abs(blx)'*abs(bly)).*sqrt(1 + 1/2*dbx.shoulder);%.*sqrt(1 + 1/2*dbx.shoulder);
                   end
                case 'coh'
                   coh(:,:,i,t,k,permi) = diag(sum(abs(blx).^2).^-.5)*csp(:,:,i,t,k)*diag(sum(abs(bly).^2).^-.5);
                    if nargout > 6
                        bias(:,:,i,t,k,permi) = sqrt(abs(blx)'.^2*abs(bly).^2)./(abs(blx)'*abs(bly)).*sqrt(1 + 1/2*dbx.shoulder);
                    end
                case {'plv','blpl','blplv','icplv'}
                    coh(:,:,i,t,k,permi) = diag(sum(abs(blx).^2).^-.5)*csp(:,:,i,t,k)*diag(sum(abs(bly).^2).^-.5);
                otherwise 
                    error('Unrecognized type %s',type)
  
            end
        end

        if nargout > 5
%             cblx = blx-repmat(mean(blx),size(blx,1),1);
%             cbly = bly-repmat(mean(bly),size(bly,1),1);
%            trf(:,:,i,t,k) = diag(sum(abs(blx).^2))\(blx'*bly);
           trf(:,:,i,t,k,permi) = (blx'*bly)/diag(sum(abs(bly).^2));
        end
            
       end
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

%%%%
function q = itercent(q)

N = 10;

for k = 1:N, 
    q = q-repmat(mean(q),[size(q,1) 1 1]); 
    q=q./abs(q); 
%    Q(:,k) = q; 
end


