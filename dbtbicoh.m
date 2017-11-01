function [out,dbx,indices] = dbtbicoh(x1,x2,varargin)

%
% [out,dbx] = dbtbicoh(x1,x2,fs,bw,varargin)
%
% Computes bicoherence with DBT.
%
% This script is under development and is not stable. Do not use it unless
% you know what you're doing.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2016

warning('This script is under development and is not stable. Do not use it unless you know what you''re doing')


dbx1 = [];
dbx2 = [];
if isa(x1,'dbt')
    if isscalar(x1)
        dbx1 = x1;        
        bw1 = dbx1.bandwidth;
        fs = dbx1.fullFS;
    else
        dbx1=x1(1);
        dbx2 =x1(2);
    end
end
    
if isscalar(x2)
    if isa(x2,'dbt')
        dbx2 = x2;
        bw2 = dbx2.bandwidth;
    elseif isnumeric(x2)
        fs = x2;
        bw1=varargin{1};
        varargin(1)=[];
        x2 = x1;
    end
elseif isnumeric(x2)
    fs = varargin{1};
    bw1 = varargin{2};
    varargin(1:2)=[];
end
    


opts = struct(...
'bw1',bw1,...
'bw2x',1,...
'symmetric',false,...
'upsampfx',0,...
'rotation','ww',...
'maxfreq',Inf,...
'type','bbb',...
'do_decomp',false,...
'decomp','svd',...
'tdopts',{{}},...
'tdrank',Inf,...
'w1lim',Inf,...
'w2lim',Inf,...
'loopchan',true,...
'normalization','awplvbc',...
'wrap',false);

if length(varargin)> 0 && isstruct(varargin{1})
    newopts = varargin{1};
    fn = fieldnames(newopts);
    for k = 1:length(fn)
        opts.(fn{k}) = newopts.(fn{k});
    end
    varargin(1) = [];
end


i = 1;     
while i < length(varargin)
  switch lower(varargin{i})

      case {'type'}
          opts.type = varargin{i+1};
          if strcmpi(opts.type,'bbb')||strcmpi(opts.type,'nnn')
              opts.bw2x = 1;
          end
          i = i+1;
      case {'symmetric'}
          opts.symmetric = varargin{i+1};
          i = i+1;
      case {'w1lim'} % Upper limit on frequency 1
          opts.w1lim= varargin{i+1};
          i = i+1;
      case {'w2lim'}
          opts.w2lim= varargin{i+1};
          i = i+1;
      case {'maxfreq','lowpass'}
          opts.maxfreq = varargin{i+1};
          i = i+1;
      case {'upsampfx'}
          opts.upsampfx = varargin{i+1};
          i = i+1;
      case {'bw2x'} % bw2 relative to bw when bandwidths differ
         opts.bw2x = varargin{i+1};
         opts.bw2=opts.bw1*opts.bw2x;
          i = i+1;
      case {'bw2'} 
          opts.bw2=varargin{i+1};
         opts.bw2x = opts.bw2/opts.bw1;
          i = i+1;
      case {'bw1'} 
         opts.bw1 = varargin{i+1};
         opts.bw2x=opts.bw2/opts.bw1;
          i = i+1;
             
    case {'svd','cp_als','cp_fastals','tucker_als'} % tensor or svd decomposition
         opts.do_decomp = varargin{i+1};
         opts.decomp = varargin{i};
          i = i+1;
 
    case {'decomp'} % tensor decomposition specified
         opts.decomp = varargin{i+1};
         opts.do_decomp = ischar(opts.decomp) && ~isempty(opts.decomp) && ~strcmpi(opts.decomp,'none');
          i = i+1;
      case {'tdrank','npc','ncomp','pckeep'} % rank of decomposition (number of components to keep for svd)
         opts.tdrank = varargin{i+1};
          i = i+1;
      case {'loopchan'} % copmute cross-bicoherece in loop rather than vectorized (for the sake of memory)
         opts.loopchan = varargin{i+1};
          i = i+1;
      case {'wrap'} % wrap the spectrum so that values are sampled modulo nyquist.
         opts.wrap = varargin{i+1};
          i = i+1;
        case {'normalization'} % copmute cross-bicoherece in loop rather than vectorized (for the sake of memory)
         opts.normalization = varargin{i+1};
          i = i+1;
            
      otherwise
        error('Unrecognized keyword %s',varargin{i})  
  end
  i=i+1;
end

opts.bw2 = opts.bw2x*opts.bw1;      
lpf = min(fs/2,opts.maxfreq*2);
nch1 = size(x1,2);
nch2 = size(x2,2);
if opts.loopchan
    
    nloop = nch2;
    nch2 = 1;
else
    nloop = 1;
%     nch2 = nch;
end

if max(nch1,nch2) > 1 && ~strcmpi(opts.decomp,'svd') && ~strcmpi(opts.decomp,'none');
   q = which(opts.decomp);
       
   while isempty(q)
        warning('Unable to find the decomposition method %s',opts.decomp)
        switch opts.decomp
            case {'cp_als','tucker_als'}
                fprintf('\nThis method requires the tensor toolbox from here:\n')
                fprintf('http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html\n')
            case {'cp_fastals'}
                fprintf('\nThis method requires the tensor toolbox from here:\n')
                fprintf('\thttp://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html\n')
                fprintf('and also the extension here: \n')
                fprintf('\thttp://www.bsp.brain.riken.jp/~phan/#tensorbox\n')               
        end
        fprintf('Do any necessary installation and add the folder to the path');
        fn = uigetdir(sprintf('Find the directory containting %s.m',opts.decomp));
        if isnumeric(fn)
            fprintf('Aborting.')
            return
        end
        addpath(fn)
        q = which(opts.decomp);
   end
end

switch lower(opts.type)
    case {'single','bbb','nnn','nnb','nbb','bnb','mi'} 
        %% The standard approach
                
                
            
                 switch lower(opts.type)
                     case {'nnb','nbb','bnb','mi'}
                         upsamptx = opts.bw2/opts.bw1-1;
                       %  upsamptx = 0;
                         upsampf2 = (opts.upsampfx+1)*(upsamptx+1)-1;
                         
                         if isscalar(x2)
                             x2 = x1;
                         end
                          dbx2 = dbt(x2,fs,opts.bw2,'remodphase',false,'upsampleFx',upsampf2,'lowpass',lpf,'gpu',false);
                     %     dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',opts.upsampfx);
                          dbx1 = dbt(x1,fs,opts.bw1,'upsampleFx',opts.upsampfx,'remodphase',false,'upsampleTx',upsamptx,'lowpass',lpf,'gpu',false);
                    
                     otherwise
                         
                         upsamptx = 1;
                         if isa(x1,'dbt')
                             dbx1 = x1;
                         else
                              dbx1 = dbt(x1,fs,opts.bw1,'upsampleFx',opts.upsampfx,'remodphase',false,'upsampleTx',upsamptx*opts.bw2x-1,'lowpass',lpf,'gpu',false);
                         end
                          if isscalar(x2)
                            dbx2 = dbx1;
                         else
                            dbx2 = dbt(x2,fs,opts.bw2,'upsampleFx',opts.upsampfx,'remodphase',false,'upsampleTx',upsamptx-1,'lowpass',lpf,'gpu',false);
                         end
                 end
                 %tol = dbx1.fullFS./dbx1.fullN/2;

         %        fstep1 = diff(dbx1.frequency(1:2));
                 fstep2 = diff(dbx2.frequency(1:2));
                 fstep3 = fstep2;
                 keept1 = dbx1.time <=dbx2.time(end) ;
                 keept2 = dbx2.time <=dbx1.time(end);%& abs(mod(dbx2.time,1./dbx1.sampling_rate))<tol/dbx1.sampling_rate;
                 keept1(sum(keept2)+1:end)=[];
                 getfreq = min(opts.maxfreq,dbx1.fullFS/2);
                
                 % Frequencies defining the region of the bispectrum to
                 % compute
                 getf = dbx1.frequency<=getfreq & dbx1.frequency>=-1e-9 & dbx1.frequency<=opts.w1lim*2; %& abs(mod(dbx1.frequency,fstep2))<tol;
                 getf2 = dbx2.frequency<=getfreq & dbx2.frequency>=-1e-9 & dbx2.frequency<=opts.w2lim;
                 % Frequency bands needed in computing the bispecturm 
                 keepf1=dbx1.frequency>=-1e-9;
                 keepf2=dbx2.frequency>=-1e-9;
                  
                 w1 = dbx1.frequency(getf);
                 w2 = dbx2.frequency(getf2);
                 if ~opts.symmetric
                      w1keep = (1:sum(getf)/2);
                      w1 = w1(w1keep);

                 else
                       w1keep = (1:sum(getf));
                 end
                 switch opts.rotation
                     case 'ww'
                         
                         %[W1,W2] = ndgrid(1:sum(getf),1:sum(getf2));   
                         if max(nch1,nch2) > 1
                            [W1,W2,WCh1,WCh2] = ndgrid(w1keep,1:opts.bw2x:sum(getf2),1:nch1,1:nch2);      
                         else
                             [W1,W2] = ndgrid(w1keep,1:opts.bw2x:sum(getf2));
                             WCh1 = 1;
                             WCh2 = 1;
                         end
                         F3 = w1(W1)+w2(W2)-w2(1);
                         
                         w2 = w2(1:opts.bw2x:sum(getf2));
%                      case 'vv'
%                          [V1,V2] = ndgrid(1:sum(getf),1:sum(getf2));        
%                         W1=(V1-V2)+1;
%                         W2=(V1+V2)-1;
%                         W3 = 2*V1-1;

                 end
%                  if opts.wrap
%                     F3 = mod(F3+dbx1.fullFS/2,dbx1.fullFS)-dbx1.fullFS/2;
%                  end
%                  W3 = (round(abs(F3)./fstep3)+1).*sign(F3);
                 W3 = round(F3./fstep3)+1;
                 %  inds = W3<=length(dbx.frequency);
                 if ~opts.wrap
                     inds = W1>0 & W2 <=sum(getf2)& abs(W3)<=sum(keepf2);
                 else
                     inds= W1>0 & W2 <=sum(getf2);
                 end
                 if ~opts.symmetric
                    inds = inds & W1 <=W2; 
                 end
                blrep = dbx1.blrep(keept1,keepf1,:);
                if nch1 > 1
                     blrep = reshape(blrep,size(blrep,1),size(blrep,2)*nch1);
                end 
%                  bspect = nan([size(W1,1) size(W1,2) nch^2]);
                 bspect = nan(size(W1));
                 
                 Ich1 =( WCh1 - 1)*sum(keepf1);
                 Ich2 =( WCh2 - 1)*sum(keepf2);
                 I1=W1 + Ich1;
                 I2=W2 + Ich2;
                 I3=W3 + Ich2;
                     if strcmpi(opts.type,'bnb')
                          I1x = I2;
                          I2 = I1;
                          I1 = I1x;
                          blrep2 = blrep;
                     elseif strcmpi(opts.type,'ppc')
                        I3 = I2; 
                        W3 = W2;
                     end
                      
                % I1=V1;I2=V2;I3=V3;
%                   if strcmpi(opts.type,'nbb') ||strcmpi(opts.type,'bbb') || strcmpi(opts.type,'nnn')
                 for chi = 1:nloop 
                     if strcmpi(opts.type,'bnb')
                        blrep = dbx2.blrep(keept2,keepf1,chi+(1:nch2)-1);
                        cblrep = conj(blrep);
                    else
                          blrep2 = dbx2.blrep(keept2,keepf2,chi+(1:nch2)-1);
                        cblrep = conj(blrep2);

                     end
                     if nch2 >1
%                          blrep = reshape(blrep,size(blrep,1),size(blrep,2)*nch);
                         blrep2 = reshape(blrep2,size(blrep2,1),size(blrep2,2)*nch2);
                         cblrep = reshape(cblrep,size(cblrep,1),size(cblrep,2)*nch2);
                     end
                     if opts.wrap
                         nyqf = min(dbx2.frequency(end),dbx2.fullFS/2);
                         wrapind = F3 >nyqf;
                         wrI = I3(wrapind);
                         [unq,unqi] = unique(F3(wrapind));
                         wrneg = (2*nyqf-unq)/fstep3+1;
                         cblrep(:,wrI(unqi))=conj(cblrep(:,wrneg));
                     end
    %                   elseif  strcmpi(opts.type,'bnb')
    %                           cblrep = conj(blrep);
    %                           blrep2 = dbx2.blrep(keept2,getf2,:);
    %                   elseif  strcmpi(opts.type,'bbb') || strcmpi(opts.type,'nnn')
    %                           cblrep = conj(blrep);
    %                           blrep2 = blrep;
    %                   end

                     bspect(inds) = sum(blrep(:,I1(inds)).*blrep2(:,I2(inds)).* cblrep(:,I3(inds)));
                    
                     NRM = nan(size(bspect));

                     bias =nan;
                     switch opts.normalization
                         case {'bicoh','standard'}
                             NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds))).^2).* sum(abs(cblrep(:,I3(inds))).^2));
                         case 'awplv'
                             NRM(inds) = sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds)).*cblrep(:,I3(inds))));
                         case 'awplvbc'
                             NRM(inds) = sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds)).*cblrep(:,I3(inds))));
                             bias = zeros(size(NRM));
                             bias(inds) = sqrt(sum((abs(blrep(:,I1(inds)).*blrep2(:,I2(inds)).*cblrep(:,I3(inds)))).^2))./sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds)).*cblrep(:,I3(inds))));
                         case 'pac'
                                   NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds))).^2).*sum(abs(blrep2(:,I2(inds))).^2.*abs(cblrep(:,I3(inds))).^2));
                     otherwise
                             error('Unrecognized normalization, %s',opts.normalization)
                     end

                     if opts.loopchan
                         bspects(:,:,:,chi) = bspect;
                         NRMs(:,:,:,chi) = NRM;
                         if ~isscalar(bias) || isnan(bias)
                             biases(:,:,:,chi) = bias;
                         end
                     end
                 end
%                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds),:)).^2).* sum(abs(cblrep(:,I3(inds),:).*blrep2(:,I2(inds),:)).^2));
             
end
if opts.loopchan
    bspect = bspects;
    NRM = NRMs;
    bias = biases;
end

out.BICOH = bspect./NRM;
out.bspect = bspect;
out.NRM = NRM;
out.bw1=opts.bw1;
out.bw2=opts.bw2;
out.bias = bias;
out.w1 = w1;
out.w2 = w2;
out.opts =opts;
out.BICOH(isnan(out.BICOH))=0;
dbx = [dbx1,dbx2];

if opts.do_decomp %%% Apply tensor or svd decomposition to the inputs
   td.Act = 0; 
  for chi = 1:nloop 
         if strcmpi(opts.type,'bnb')
            blrep = dbx2.blrep(keept2,keepf1,chi+(1:nch2)-1);
            cblrep = conj(blrep);
        else
              blrep2 = dbx2.blrep(keept2,keepf2,chi+(1:nch2)-1);
            cblrep = conj(blrep2);

         end

        if max(nch1,nch2) ==1 || strcmpi(opts.decomp,'svd') 
            if max(nch1,nch2) >1
                BICOH = reshape(out.BICOH,size(out.BICOH(:,:,1)).*[1 nch1*nch2]);
            else
                BICOH = out.BICOH;
            end
            [u,l,v] = svd(BICOH); 
            getn = min(sum(diag(l)>0),opts.tdrank);
            [uind,vind] = ndgrid(1:size(BICOH,1),1:size(BICOH,2),1:nch1,1:nch2);
            U = u(uind(inds),1:getn);
            V = v(vind(inds),1:getn);

    %     switch opts.type
    %         case 'nbb'

    %            A =  (blrep(:,I1(inds)).*blrep2(:,I2(inds)).* cblrep(:,I3(inds)))*(repmat(NRM(inds).^-1,1,getn).*(V.*conj(U))*diag(diag(l).^-.5));
    %         otherwise
    %             warning('SVD not yet implemented for case %s',opts.type)
    %             A = nan;
    %     end
         %   tdcomp.Act = A;
            tdcomp.time = dbx1.time;
            tdcomp.u = u(:,1:getn);
            tdcomp.l = diag(l);
            if max(nch1,nch2) ==1
                tdcomp.v = v(:,1:getn);
            else
               tdcomp.v = reshape(v(:,1:getn), [size(v,1)/(nch1*nch2)  nch1 nch2 getn]);
            end
            UMAT= (V.*conj(U))*diag(diag(l(1:getn,1:getn)).^-.5);
        else
            %%% Apply tensor decomposition
            getn = opts.tdrank;
            decomp = str2func(opts.decomp); 
            tdcomp = decomp(tensor(out.BICOH),opts.tdrank,opts.tdopts{:});
            [u1ind,u2ind,ch1ind,ch2ind] = ndgrid(1:size(out.BICOH,1),1:size(out.BICOH,2),1:nch1,chi+(1:nch2)-1);
             UMAT = conj((tdcomp.U{1}(u1ind(inds),:).*tdcomp.U{2}(u2ind(inds),:).*tdcomp.U{3}(ch1ind(inds),:).*tdcomp.U{4}(ch2ind(inds),:))*diag(sqrt(tdcomp.lambda).^-1));
            tdcomp = struct('tensor',tdcomp);
        end
          % Time activation
           tdcomp.Act = td.Act + (blrep(:,I1(inds)).*blrep2(:,I2(inds)).* cblrep(:,I3(inds)))*(repmat(NRM(inds).^-1,1,getn).*UMAT);
            tdcomp.time = dbx1.time;

  end
else
         tdcomp = [];
end

out.nyqf = dbx1.fullFS/2;
out.tdcomp = tdcomp;

if nargout > 2
    indices = [I1(inds),I2(inds),I3(inds)];
end