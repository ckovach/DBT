function [BICOH,bspect,w1,w2,NRM,dbx1,dbx2,A] = dbtbicoh(x,fs,bw,varargin)

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

%if nargin < 3
    bw2x = 4;
%end
downsamp2 = true; 
symmetric = false;
upsampfx = 0;
rotation = 'ww';
maxfreq = Inf;
type = 'bbb';

  i = 1;   
  
while i < length(varargin)
  switch lower(varargin{i})

      case {'type'}
          type = varargin{i+1};
          if strcmpi(type,'bbb')||strcmpi(type,'nnn')
              bw2x = 1;
          end
          i = i+1;
      case {'symmetric'}
          symmetric = varargin{i+1};
          i = i+1;
      case {'maxfreq'}
          maxfreq = varargin{i+1};
          i = i+1;
      case {'upsampfx'}
          upsampfx = varargin{i+1};
          i = i+1;
      case {'bw2x'} % bw2 relative to bw when bandwidths differ
         bw2x = varargin{i+1};
          i = i+1;
      case {'bw2'} % bw2 relative to bw when bandwidths differ
         bw2x = varargin{i+1}/bw;
          i = i+1;
             
        case {'svd'} % bw2 relative to bw when bandwidths differ
         do_svd = varargin{i+1};
          i = i+1;
           
      otherwise
        error('Unrecognized keyword %s',varargin{i})  
  end
  i=i+1;
end

bw2 = bw2x*bw;      
lpf = min(fs/2,maxfreq*2);
nch = size(x,2);
switch lower(type)
    case {'single','bbb','nnn','nnb','nbb','bnb'} 
        %% The standard approach
                
                
            
                 switch lower(type)
                     case {'nnb','nbb','bnb'}
                         upsamptx = bw2/bw-1;
                       %  upsamptx = 0;
                         upsampf2 = (upsampfx+1)*(upsamptx+1)-1;
                          dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',upsampf2,'lowpass',lpf);
                     %     dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',upsampfx);
                          dbx1 = dbt(x,fs,bw,'upsampleFx',upsampfx,'remodphase',false,'upsampleTx',upsamptx,'lowpass',lpf);
                    
                     otherwise
                         
                         upsamptx = 1;
                          dbx1 = dbt(x,fs,bw,'upsampleFx',upsampfx,'remodphase',false,'upsampleTx',upsamptx-1,'lowpass',lpf);
                          dbx2 = dbx1;
                 end
                 tol = dbx1.fullFS./dbx1.fullN/2;

         %        fstep1 = diff(dbx1.frequency(1:2));
                 fstep2 = diff(dbx2.frequency(1:2));
                 fstep3 = fstep2;
                 keept1 = dbx1.time <=dbx2.time(end) ;
                 keept2 = dbx2.time <=dbx1.time(end);%& abs(mod(dbx2.time,1./dbx1.sampling_rate))<tol/dbx1.sampling_rate;
                 keept1(sum(keept2)+1:end)=[];
                 getfreq = min(maxfreq,dbx1.fullFS/2);
                
                 getf = dbx1.frequency<=getfreq & dbx1.frequency>=-1e-9; %& abs(mod(dbx1.frequency,fstep2))<tol;
                 getf2 = dbx2.frequency<=getfreq & dbx2.frequency>=-1e-9;
                 w1 = dbx1.frequency(getf);
                 w2 = dbx2.frequency(getf2);
                 if ~symmetric
                      w1keep = (1:sum(getf)/2);
                      w1 = w1(w1keep);

                 else
                       w1keep = (1:sum(getf));
                 end
                 switch rotation
                     case 'ww'
                         
                         %[W1,W2] = ndgrid(1:sum(getf),1:sum(getf2));   
                         if nch > 1
                            [W1,W2,WCh1,WCh2] = ndgrid(w1keep,1:bw2x:sum(getf2),1:nch,1:nch);      
                         else
                             [W1,W2] = ndgrid(w1keep,1:bw2x:sum(getf2));
                             WCh1 = 1;
                             WCh2 = 1;
                         end
                         W3 = round((w1(W1)+w2(W2)-w2(1))./fstep3)+1;
                         
                         w2 = w2(1:bw2x:sum(getf2));
                     case 'vv'
                         [V1,V2] = ndgrid(1:sum(getf),1:sum(getf2));        
                        W1=(V1-V2)+1;
                        W2=(V1+V2)-1;
                        W3 = 2*V1-1;

                 end
             
                 %  inds = W3<=length(dbx.frequency);
                 inds = W1>0 & W2 <=sum(getf2)& W3<=sum(getf2);
                 if ~symmetric
                    inds = inds & W1 <=W2; 
                 end
                 blrep = dbx1.blrep(keept1,getf,:);
%                  bspect = nan([size(W1,1) size(W1,2) nch^2]);
                 bspect = nan(size(W1));
                 
                 Ich1 =( WCh1 - 1)*sum(getf);
                 Ich2 =( WCh2 - 1)*sum(getf2);
                 I1=W1 + Ich1;
                 I2=W2 + Ich2;
                 I3=W3 + Ich2;
                      
                % I1=V1;I2=V2;I3=V3;
%                   if strcmpi(type,'nbb') ||strcmpi(type,'bbb') || strcmpi(type,'nnn')
                  
                 if strcmpi(type,'bnb')
                      W1x = W2;
                      W2 = W1;
                      W1 = W1x;
                      blrep2 = blrep;
                    blrep = dbx2.blrep(keept2,getf2,:);
                    cblrep = conj(blrep);
                else
                      blrep2 = dbx2.blrep(keept2,getf2,:);
                    cblrep = conj(blrep2);
                
                 end
                 if nch >1
                     blrep = reshape(blrep,size(blrep,1),size(blrep,2)*nch);
                     blrep2 = reshape(blrep2,size(blrep2,1),size(blrep2,2)*nch);
                     cblrep = reshape(cblrep,size(cblrep,1),size(cblrep,2)*nch);
                 end
%                   elseif  strcmpi(type,'bnb')
%                           cblrep = conj(blrep);
%                           blrep2 = dbx2.blrep(keept2,getf2,:);
%                   elseif  strcmpi(type,'bbb') || strcmpi(type,'nnn')
%                           cblrep = conj(blrep);
%                           blrep2 = blrep;
%                   end

                 bspect(inds) = sum(blrep(:,I1(inds)).*blrep2(:,I2(inds)).* cblrep(:,I3(inds)));
                    
                 NRM = nan(size(bspect));
                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds)).*blrep2(:,I2(inds))).^2).* sum(abs(cblrep(:,I3(inds))).^2));
%                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds),:)).^2).* sum(abs(cblrep(:,I3(inds),:).*blrep2(:,I2(inds),:)).^2));
             
end
   BICOH = bspect./NRM;
                     
   
   BICOH(isnan(BICOH))=0;
   
 if do_svd 
    [u,l,v] = svd(BICOH); 
    getn = sum(diag(l)>0);
    [uind,vind] = ndgrid(1:size(BICOH,1),1:size(BICOH,2));
    U = u(uind(inds),1:getn);
    V = v(vind(inds),1:getn);
    switch type
        case 'nbb'
                
            A =  (blrep(:,I1(inds),:).*blrep2(:,I2(inds),:).* cblrep(:,I3(inds),:))*(repmat(NRM(inds).^-1,1,getn).*(V.*conj(U))*diag(diag(l).^-.5));
        otherwise
            warning('SVD not yet implemented for case %s',type)
            A = nan;
 end
 else
     A = nan;
 end