function [out,dbx] = dbtbicoh(x,fs,bw,varargin)

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

opts = struct(...
'bw',bw,...
'bw2x',4,...
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
'w2lim',Inf);

if nargin > 3 && isstruct(varargin{1})
    newopts = varargin{1};
    fn = fieldnames(newopts);
    for k = 1:lengt(fn)
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
      case {'maxfreq'}
          opts.maxfreq = varargin{i+1};
          i = i+1;
      case {'upsampfx'}
          opts.upsampfx = varargin{i+1};
          i = i+1;
      case {'bw2x'} % bw2 relative to bw when bandwidths differ
         opts.bw2x = varargin{i+1};
          i = i+1;
      case {'bw2'} 
         opts.bw2x = varargin{i+1}/opts.bw;
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
            
      otherwise
        error('Unrecognized keyword %s',varargin{i})  
  end
  i=i+1;
end

bw2 = opts.bw2x*opts.bw;      
lpf = min(fs/2,opts.maxfreq*2);
nch = size(x,2);

if nch > 1 && ~strcmpi(opts.decomp,'svd') && ~strcmpi(opts.decomp,'none');
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
    case {'single','bbb','nnn','nnb','nbb','bnb'} 
        %% The standard approach
                
                
            
                 switch lower(opts.type)
                     case {'nnb','nbb','bnb'}
                         upsamptx = bw2/opts.bw-1;
                       %  upsamptx = 0;
                         upsampf2 = (opts.upsampfx+1)*(upsamptx+1)-1;
                          dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',upsampf2,'lowpass',lpf);
                     %     dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',opts.upsampfx);
                          dbx1 = dbt(x,fs,opts.bw,'upsampleFx',opts.upsampfx,'remodphase',false,'upsampleTx',upsamptx,'lowpass',lpf);
                    
                     otherwise
                         
                         upsamptx = 1;
                          dbx1 = dbt(x,fs,opts.bw,'upsampleFx',opts.upsampfx,'remodphase',false,'upsampleTx',upsamptx-1,'lowpass',lpf);
                          dbx2 = dbx1;
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
                         if nch > 1
                            [W1,W2,WCh1,WCh2] = ndgrid(w1keep,1:opts.bw2x:sum(getf2),1:nch,1:nch);      
                         else
                             [W1,W2] = ndgrid(w1keep,1:opts.bw2x:sum(getf2));
                             WCh1 = 1;
                             WCh2 = 1;
                         end
                         W3 = round((w1(W1)+w2(W2)-w2(1))./fstep3)+1;
                         
                         w2 = w2(1:opts.bw2x:sum(getf2));
                     case 'vv'
                         [V1,V2] = ndgrid(1:sum(getf),1:sum(getf2));        
                        W1=(V1-V2)+1;
                        W2=(V1+V2)-1;
                        W3 = 2*V1-1;

                 end
             
                 %  inds = W3<=length(dbx.frequency);
                 inds = W1>0 & W2 <=sum(getf2)& W3*fstep3<=opts.maxfreq;
                 if ~opts.symmetric
                    inds = inds & W1 <=W2; 
                 end
                 blrep = dbx1.blrep(keept1,keepf1,:);
%                  bspect = nan([size(W1,1) size(W1,2) nch^2]);
                 bspect = nan(size(W1));
                 
                 Ich1 =( WCh1 - 1)*sum(keepf1);
                 Ich2 =( WCh2 - 1)*sum(keepf2);
                 I1=W1 + Ich1;
                 I2=W2 + Ich2;
                 I3=W3 + Ich2;
                      
                % I1=V1;I2=V2;I3=V3;
%                   if strcmpi(opts.type,'nbb') ||strcmpi(opts.type,'bbb') || strcmpi(opts.type,'nnn')
                  
                 if strcmpi(opts.type,'bnb')
                      I1x = I2;
                      I2 = I1;
                      I1 = I1x;
                      blrep2 = blrep;
                    blrep = dbx2.blrep(keept2,keepf1,:);
                    cblrep = conj(blrep);
                else
                      blrep2 = dbx2.blrep(keept2,keepf2,:);
                    cblrep = conj(blrep2);
                
                 end
                 if nch >1
                     blrep = reshape(blrep,size(blrep,1),size(blrep,2)*nch);
                     blrep2 = reshape(blrep2,size(blrep2,1),size(blrep2,2)*nch);
                     cblrep = reshape(cblrep,size(cblrep,1),size(cblrep,2)*nch);
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
                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds))).^2).* sum(abs(cblrep(:,I3(inds))).^2));
%                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds),:)).^2).* sum(abs(cblrep(:,I3(inds),:).*blrep2(:,I2(inds),:)).^2));
             
end
out.BICOH = bspect./NRM;
out.bspect = bspect;
out.NRM = NRM;
out.w1 = w1;
out.w2 = w2;
out.opts =opts;
out.BICOH(isnan(out.BICOH))=0;
dbx = [dbx1,dbx2];

if opts.do_decomp %%% Apply tensor or svd decomposition to the inputs
    if nch ==1 || strcmpi(opts.decomp,'svd') 
        if nch >1
            BICOH = reshape(out.BICOH,size(out.BICOH(:,:,1)).*[1 nch^2]);
        else
            BICOH = out.BICOH;
        end
        [u,l,v] = svd(BICOH); 
        getn = min(sum(diag(l)>0),opts.tdrank);
        [uind,vind] = ndgrid(1:size(BICOH,1),1:size(BICOH,2),1:nch,1:nch);
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
        if nch ==1
            tdcomp.v = v(:,1:getn);
        else
           tdcomp.v = reshape(v(:,1:getn), [size(v,1)/nch^2  nch nch getn]);
        end
        UMAT= (V.*conj(U))*diag(diag(l(1:getn,1:getn)).^-.5);
    else
        %%% Apply tensor decomposition
        getn = opts.tdrank;
        decomp = str2func(opts.decomp); 
        tdcomp = decomp(tensor(out.BICOH),opts.tdrank,opts.tdopts{:});
        [u1ind,u2ind,ch1ind,ch2ind] = ndgrid(1:size(out.BICOH,1),1:size(out.BICOH,2),1:nch,1:nch);
         UMAT = conj((tdcomp.U{1}(u1ind(inds),:).*tdcomp.U{2}(u2ind(inds),:).*tdcomp.U{3}(ch1ind(inds),:).*tdcomp.U{4}(ch2ind(inds),:))*diag(sqrt(tdcomp.lambda).^-1));
        tdcomp = struct('tensor',tdcomp);
    end
      % Time activation
       tdcomp.Act =  (blrep(:,I1(inds)).*blrep2(:,I2(inds)).* cblrep(:,I3(inds)))*(repmat(NRM(inds).^-1,1,getn).*UMAT);
        tdcomp.time = dbx1.time;
 
else
     tdcomp = [];
end
out.tdcomp = tdcomp;

