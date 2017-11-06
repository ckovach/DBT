function [out,pspindices] = pspect2(S,varargin)

%[psp,pspindices] = pspect2({S1,S2,..,Sk},{f1,f2,...,fk},order,[options)
%
% General function for computing polyspectra and cross polyspectra of the 
% given  order from a cell array of time x frequency matrices.
%
% Input arguments:
%       S:  1 x order cell array with matrices of time-frequency
%           coefficients as time x frequency. If the length of S is less than
%           order, then the last element of S is reduplicated to extend S to
%           length order.
%       f:  1 x order (or length(S) if it is less than order) cell array of 
%            vectors with frequency labels for the columns of S.
%          
%   order:  Order of the polyspectrum.
%
%
%Options: Options may either be specified as pairs of keywords and values, 
%         i.e. pspect(...,'keyword',value,...) or as one or more structs
%         with  struct.keyword = value.
%   Avalailable options are:      
%
%   lowpass: limit the range of frequencies for each axis to values less-
%           than or equal to this. This may be specified as a scalar or as
%           a vector of order-1 length, which applies a separate limit for
%           each dimension.
%   highpass: limit to frequencies above this value, etc. 
%   maxfreq: Sum of frequencies across all dimensions must be less-than-or
%            -equal to this value.
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
%   Input may alternatively be a 1 x order dbt object array. Note that each
%   has to have the same sampling rate, meaning that if bandwidths differ,
%   upsampleTx has to be adjusted accordingly. 
%
%
% See also DBT, PSPECT

% C. Kovach 2017

options.lowpass= Inf;
options.maxfreq= Inf;
options.highpass= 0;
options.normalization = 'awplv';
options.full_range = false; % Add negative frequencies if they are not already included
options.min_range = false; % Return only one symmetry region
options.symmetrize = false;
options.round_freq = true; % Round to the nearest frequency band if necessary.
options.tolerance = []; % Rounding tolerance (defaults to min(diff(f))).
%options.real_signal=true;

axes_interchangeable = isnumeric(S) || length(S)==1;

if isa(S,'dbt')
    
    order = varargin{1};
    dbx = S;
    dbx(end+1:order)=dbx(end);
    
    fs = {dbx.frequency};
    mint = min(arrayfun(@(d)d.time(end),dbx));
    S={};
    for k = 1:order
        S{k} = dbx(k).blrep(dbx(k).time<=mint,:);
    end
    varargin(1)=[];
else
    f = varargin{1};
    order=varargin{2};
    varargin(1:2)=[];
    if ~iscell(S)
        S = repmat({S},1,order);
    else
        S(end+1:order)=S(end);
    end
    if ~iscell(f)
         fs = repmat({f},1,order-1);
    else
        fs=f;
        fs(end+1:order)=f;
    end

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

options.full_range = options.full_range && ~options.min_range;

if isscalar(options.lowpass)
    options.lowpass = options.lowpass*ones(1,order);
end

if isscalar(options.highpass)
    options.highpass = options.highpass*ones(1,order);
end



if isempty(options.tolerance)
    options.tolerance = cellfun(@(f)min(diff(sort(f))),fs); 
end

fs = arrayfun(@(f,highpass,lowpass)f{1}(abs(f{1})>=highpass & abs(f{1})<=lowpass),fs,options.highpass,options.lowpass,'uniformoutput',false);


for k =1:order
    Sk=S{k};
    if options.full_range

        sindx = 1:size(Sk,2);
        f=fs{k};
        [fs{k},srtf] = unique([f,-f]);
        sifull = [sindx,sindx];
        sconj{k} = [false(size(sindx)),true(size(sindx))];
        sconj{k}=sconj{k}(srtf);

        resortindex{k} = sifull(srtf);
        Sk = Sk(:,resortindex{k});
        Sk(:,sconj{k}) = conj(Sk(:,sconj{k}));
        S{k}=Sk;
    else
        resortindex{k} = 1:size(Sk,2); %#ok<*AGROW>
        sconj{k} = false(size(resortindex{k}));
    end
end

W = fs(1:order-1);
[W{:}] = ndgrid(fs{1:order-1});
% W{order} = -sum(cat(order,W{:}),order);

WW = cellfun(@(x)x(:),W,'uniformoutput',false);
WW = [WW{:}];

if options.symmetrize
    n1 = 1/(order-1);
    main_dim = 1;
    other_dims = setdiff(1:order-1,main_dim);
    WW(:,other_dims) = WW(:,other_dims)-n1*repmat(WW(:,main_dim),1,size(WW,2)-1);
end

WW(:,order) = -sum(WW,2);

if axes_interchangeable
    WW = sort(WW,2);
end

% unique combinations only
[wunq] = unique(WW,'rows');
%wunq = WW(wunqi,:);
mxf = cellfun(@(f)max(abs(f)),fs);
wunq(any(abs(wunq)>repmat(mxf,size(wunq,1),1),2) | abs(wunq(:,order))>options.maxfreq,:)=[];
[ism,indx] = ismember(WW,wunq,'rows');
% if options.min_range
%     [unqindx,unqi]=unique(indx);
%     indx(:)=0;
%     indx(unqi)=unqindx;
%     ism(:)=false;
%     ism(unqi)=true;
% end
    
if axes_interchangeable
    [~,cindx] = ismember(sort(-WW(~ism,:),2),wunq,'rows');
    indx(~ism)=cindx;
end



PS = 1;
NORM = 1;
for k = 1:order
    [fism,findx] = ismember(round(wunq(:,k)./options.tolerance(k)),round(fs{k}./options.tolerance(k)));

    [~,findx(~fism)] = ismember(round(-wunq(~fism,k)./options.tolerance(k)),round(fs{k}./options.tolerance(k)));
    
    fisms(:,k)=fism;
    findxs(:,k)=findx;
    
    F = S{k}(:,findx);
    F(:,~fism) = conj(F(:,~fism));
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

if options.min_range
    [uindx,uindxi] = unique(indx);
    indx(:)=0;
    indx(uindxi)=uindx;
end

indx(indx==0)=length(psp);

    
rmat = reshape(indx,size(W{1}));


out.pspect = psp(rmat);
out.fs = fs(1:2);

switch options.normalization
    case  {'awplv'}
        nrm = sum(NORM);
        nrm(end+1)=1;
        out.normalization = nrm(rmat);
        bias = sqrt(sum(NORM.^2)./sum(NORM).^2);
        bias(end+1)=0;
        out.bias = bias(rmat);
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

out.options = options;

[unq,unqi] = unique(rmat(:));
lin(unq) = unqi;
out.lin = lin(:);
out.square = rmat;


if nargout >1
    if axes_interchangeable
        pspindices.findex = resortindex{1}(findxs);
        pspindices.conjugate = xor(fisms,sconj{1}(findxs));
    else
        for k =1:order
            pspindices.findex(:,k) = resortindex{k}(findxs(:,k));
            pspindices.conjugate = xor(fisms(:,k),sconj{k}(findxs(:,k))');
        end
    end
    pspindices.reconmat = rmat;
        
end

