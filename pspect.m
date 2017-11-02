function [out,pspindices] = pspect(S,varargin)

%[psp,pspindices] = pspect(S,f,order,[options])
%
% General function for computing polyspectra of the given  order from an
% input time x frequency matrix.
%
% Input arguments:
%       S:  Matrix of time-frequency coefficients as time x frequency
%       f:  Frequency labels for the columns of S.
%   order:  Order of the polyspectrum.
%
%
%Options: Options may either be specified as pairs of keywords and values, 
%         i.e. pspect(...,'keyword',value,...) or as one or more structs
%         with  struct.keyword = value.
%   Avalailable options are:      
%
%   lowpass: limit the range of frequencies for each axis to values less-
%           than or equal to this. This maybe specified as a scalar or as
%           a vector of order-1 length, which applies a separate limit for
%           each dimension.
%   highpass: limit to frequencies above this value, etc. 
%   maxfreq: Sum of frequencies across all dimensions must be less-than-or
%            -equal to this value.
%   full_range: Add negative frequencies if they are not already included
%
%
% Output arguments:
%       psp: Struct with the following fields:
%          .pspect: Unnormalized polyspectrum as an (order-1)-dimensional
%                   matrix;
%          .fs:     Frequuency labels for the dimensions of pspect
%          .options Options struct.
%
%       pspindices: Struct containing indices into the original data with the
%            following fields:
%           .findex: index of the frequency for each unqiue term in the
%                    estimate.
%           .conjugate: Terms for which the complex conjugate is taken.
%           .reconmat: reconstruct into the same shape and size as psp.pspect. 
%
%
% [psp,pspindices] = pspect(dbx,order,[options)
% 
%   Input may alternatively be a dbt object.
%
%
% See also DBT

% C. Kovach 2017

options.lowpass= Inf;
options.maxfreq= Inf;
options.highpass= 0;
options.normalization = 'awplv';
options.full_range = true; % Add negative frequencies if they are not already included
options.symmetrize = false;
options.round_freq = true; % Round to the nearest frequency band if necessary.
options.tolerance = []; % Rounding tolerance (defaults to min(diff(f))).
options.principal_domain = false; %Only return values in the principal domain.
options.getbias = true;
%options.real_signal=true;


if isa(S,'dbt')
    f = S.frequency;
    S = S.blrep;
    order = varargin{1};
    varargin(1)=[];
else
    f = varargin{1};
    order=varargin{2};
    varargin(1:2)=[];
end
if isempty(options.tolerance)
    options.tolerance = min(diff(sort(f))); 
end


optfld = fieldnames(options);
i = 1;
while i <length(varargin)
    if isstruct(varargin{i})
        fldn = fieldnames(varargin{i});
        for k = 1:length(fldn)
            if ~ismember(fldn{k},optfld)
                error('Unrecognized option, %s',fldn{k})
            end
            options.(fldn{k}) = varargin{i}.(fldn{k});
        end
        i = i-1;
    elseif ismember(varargin{i},optfld)
        options.(varargin{i})=varargin{i+1};
    else
        error('Unrecognized option, %s',varargin{i});
    end
    i = i+2;
end

if isscalar(options.lowpass)
    options.lowpass = options.lowpass*ones(1,order-1);
end

if isscalar(options.highpass)
    options.highpass = options.highpass*ones(1,order-1);
end

if options.full_range
    sindx = 1:size(S,2);
    [f,srtf] = unique([f,-f]);
    sifull = [sindx,sindx];
    sconj = [false(size(sindx)),true(size(sindx))];
    sconj=sconj(srtf);
    
    resortindex = sifull(srtf);
    S = S(:,resortindex);
    S(:,sconj) = conj(S(:,sconj));
    
else
    resortindex = 1:size(S,2);
    sconj = false(size(resortindex));
end


% if ~options.symmetrize
     fs = repmat({f},1,order-1);
% else
%     n1 = 1/(order-1);
%     fs = {f};
%     for k = 2:order-1
%         fs{k} = f-n1*f;
%     end
% end
fs = arrayfun(@(f,highpass,lowpass)f{1}(abs(f{1})>=highpass & abs(f{1})<=lowpass),fs,options.highpass,options.lowpass,'uniformoutput',false);


W = fs;
[W{:}] = ndgrid(fs{:});
% W{order} = -sum(cat(order,W{:}),order);

principal_domain=true;
for k = 2:length(W)
    principal_domain = principal_domain & W{k-1}<=W{k};
end

if ~options.symmetrize
    WW = cellfun(@(x)x(:),W,'uniformoutput',false);
else
    n1 = 1/(order-1);
    WW = arrayfun(@(x,a)x{1}(:)-a*W{1}(:),W,n1*(1:order-1>1),'uniformoutput',false);
end

WW = [WW{:}];
WW(:,order) = -sum(WW,2);
WW = sort(WW,2);

% unique combinations only
if options.principal_domain
    [wunq] = unique(WW(principal_domain(:),:),'rows');
else
       [wunq] = unique(WW,'rows');
end
   
wunq(any(abs(wunq)>max(abs(f)),2) | abs(wunq(:,order))>options.maxfreq,:)=[];
[ism,indx] = ismember(WW,wunq,'rows');
[cpart,cindx] = ismember(sort(-WW(~ism,:),2),wunq,'rows');
indx(~ism)=cindx;

[fism,findx] = ismember(round(wunq./options.tolerance),round(f./options.tolerance));

[~,findx(~fism)] = ismember(round(-wunq(~fism)./options.tolerance),round(f./options.tolerance));

PS = 1;
NORM = 1;
for k = 1:order
    
    F = S(:,findx(:,k));
    F(:,~fism(:,k)) = conj(F(:,~fism(:,k)));
    PS = PS.*F;
    
    
    switch options.normalization
        case {'awplv'}
            NORM = NORM.*abs(F);
        case {'polycoh','coh','coherence','polycoherence','bicoherence'}
            options.normalization = 'polycoh';
            if k<order
                NORM = NORM.*F;
            else
                NORM(:,:,end+1) = abs(F).^2;
            end
        case 'rms'
            NORM = NORM.*sum(abs(F).^2);
        case {'','none'}
        otherwise
            error('Unrecognized normalization %s',options.normalization)
    end

end

psp = sum(PS);
psp(end+1)=0;
indx(indx==0)=length(psp);
if options.principal_domain
   indx(~principal_domain(:))=length(psp); 
end

rmat = reshape(indx,size(W{1}));


out.pspect = psp(rmat);
out.fs = fs;

switch options.normalization
    case  {'awplv'}
        nrm = sum(NORM);
        if options.getbias
            nrm2 = sum(NORM.*NORM);
            bias = sqrt(nrm2./nrm.^2);
             bias(end+1)=0;
            out.bias = bias(rmat);
        end
        nrm(end+1)=1;
        out.normalization = nrm(rmat);
    case 'polycoh'        
        nrm = sqrt(sum(abs(NORM(:,:,1)).^2).*sum(NORM(:,:,2)));
        nrm(end+1)=1;
        out.normalization = nrm(rmat);
        out.bias=nan;
    case 'rms'        
        nrm = sqrt(NORM);
        nrm(end+1)=1;
        out.normalization = nrm(rmat);
        out.bias = nan;
    case {'','none'}
        out.normalization = 1;
  
end

out.options = options;tic

if nargout >1
    pspindices.findex = resortindex(findx);
    pspindices.conjugate = xor(fism,sconj(findx));
    pspindices.reconmat = rmat;
    
    
    pspindices.principal_domain=principal_domain;
    
end

