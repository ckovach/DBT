function [BICOH,bspect,w,w2,NRM,dbx1,dbx2] = dbtbicoh(x,fs,bw,varargin)

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
symmetric = true;
upsampfx = 1;
rotation = 'ww';

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
      case {'bw2x'} % bw2 relative to bw when bandwidths differ
         bw2x = varargin{i+1};
          i = i+1;
             
      otherwise
        error('Unrecognized keyword %s',varargin{i})  
  end
  i=i+1;
end

if nargin < 3
    bw = 1;      
end
bw2 = bw2x*bw;      

switch lower(type)
    case 'fcov'
        %%
%        upsamptx = round(min(300,fs)/(2*bw));
%        upsamptx = 0;
                upsamptx = round(min(300,fs)/(2*bw));

        %%% Compute DBT
        dbx1 = dbt(x,fs,bw,'upsampleTx',upsamptx-1,'upsampleFx',upsampfx,'remodphase',true);
       %%% Anti-aliasing filter and decimation
        getfreq = dbx1.sampling_rate/2;
        nadj = sum(dbx1.time<=(length(x)-1)./fs);
        %%% Get lowpass and hilbert filtered data
        fxrs = fft(dbx1.signal(dbx1.frequency < getfreq,false));
        
        %%% Resample to match dbt sampling rate as nearly as possible
        fxrs(nadj+1:end)=[];
        xrs = real(ifft(fxrs));
        
        
        % decim = fs./dbx.sampling_rate;
        % xrs = x(round(decim/2):decim:end,:);
        
        %%% Weight with original signal 
        blrep = dbx1.blrep(1:length(xrs),dbx1.frequency<getfreq);
        w = dbx1.frequency(dbx1.frequency<getfreq);
        sz = size(blrep);
        sz(end+1:3)=1;
        blrep = reshape(blrep,[sz(1) sz(2)*sz(3)]);
        Wblrep = reshape(repmat(permute(xrs,[1 3 2]),1,length(w)),[sz(1) sz(2)*sz(3)]).*blrep;
        %%% Estimate bispectrum using the frequency-domain covariance
        %%% method (Kovach 2016)
        bspect = conj(blrep)'*Wblrep;
              [W1,W2] = ndgrid(1:length(w),1:length(w));
                W3 = W1+W2-1;
                inds = W3<=size(blrep,2);
         normalization ='std';
        switch normalization
%             case 'sig'
%                  %%% Signal to generate normalizing term of bicoherence
%                  dbxnrm = dbx;
%                 dbxnrm.blrep = dbxnrm.blrep;
%                 dbxnrm.remodphase=true;
%                 fnrmxrs = fft(dbxnrm.signal(dbx.frequency < getfreq,true));
%                 fnrmxrs(nadj+1:end)=[];
%                 nrmxrs = real(ifft(fnrmxrs));
%                 
%                 blrepnrm = dbxnrm.blrep;%.*dbxnrm.remodulator;
%                 blrepnrm = blrepnrm(1:length(nrmxrs),dbxnrm.frequency<getfreq);
%                 blrepnrm = reshape(blrepnrm,[sz(1) sz(2)*sz(3)]);
%                 Wblreprnm = reshape(repmat(permute(nrmxrs,[1 3 2]),1,length(w)),[sz(1) sz(2)*sz(3)]).*blrepnrm;
%                 %%% Estimate bispectrum using the frequency-domain covariance
%                 %%% method (Kovach 2016)
%                 NRM = conj(blrepnrm)'*Wblreprnm;
%               case 'cube'
%                 %%% Normalization with the root product of the sum of cubed energy
%                 %%% in each band        
%                 SC =sum(abs(blrep).^3);
%                  NRM = nan(size(bspect));
%                 NRM(inds) = (SC(W1(inds)).*SC(W2(inds)).*SC(W3(inds))).^(1/3);     
            case 'std'
                 SS = sum(abs(blrep).^2);
                 SSQ = (abs(blrep').^2*abs(blrep).^2);
                 NRM = nan(size(bspect));
                 NRM(inds) = sqrt(SSQ(inds).*SS(W3(inds))');
             case 'smstd'
                 SS = convn(sum(abs(blrep).^2),ones(1,2));
                 SSQ = (abs(blrep').^2*abs(blrep).^2);
                 NRM = nan(size(bspect));
                 NRM(inds) = sqrt(SSQ(inds).*SS(W3(inds))');
%             case 'band'
%                  dbp = dbt(abs(blrep).^2,dbx.sampling_rate,bw,'upsampleFx',upsampfx);
%                  MSS = squeeze(sum(abs(dbp.blrep(:,dbp.frequency<getfreq,:)).^2));
%                  XSS = sum(abs(blrep).^2);
%                  NRM = nan(size(bspect));
%                  NRM(inds) = sqrt(MSS(inds).*XSS(W3(inds))');
        end
%        Wnrm = reshape(repmat(permute(abs(xrs),[1 3 2]),1,length(w)),[sz(1) sz(2)*sz(3)]).*abs(blrep);
     
    case {'single','bbb','nnn','nnb','nbb','bnb'} 
        %% The standard approach
                
                
            
                 switch lower(type)
                     case {'nnb','nbb','bnb'}
                         upsamptx = bw2/bw-1;
                       %  upsamptx = 0;
                         upsampf2 = (upsampfx+1)*(upsamptx+1)-1;
                          dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',upsampf2);
                     %     dbx2 = dbt(x,fs,bw2,'remodphase',false,'upsampleFx',upsampfx);
                          dbx1 = dbt(x,fs,bw,'upsampleFx',upsampfx,'remodphase',false,'upsampleTx',upsamptx);
                    
                     otherwise
                         
                         upsamptx = 1;
                          dbx1 = dbt(x,fs,bw,'upsampleFx',upsampfx,'remodphase',false,'upsampleTx',upsamptx-1);
                          dbx2 = dbx1;
                 end
                 tol = dbx1.fullFS./dbx1.fullN/2;

         %        fstep1 = diff(dbx1.frequency(1:2));
                 fstep2 = diff(dbx2.frequency(1:2));
                 fstep3 = fstep2;
                 keept1 = dbx1.time <=dbx2.time(end) ;
                 keept2 = dbx2.time <=dbx1.time(end)& abs(mod(dbx2.time,1./dbx1.sampling_rate))<tol/dbx1.sampling_rate;
                 getfreq = dbx1.fullFS/2;
                
                 getf = dbx1.frequency<=getfreq & dbx1.frequency>=0 & abs(mod(dbx1.frequency,fstep2))<tol;
                 getf2 = dbx2.frequency<=getfreq & dbx2.frequency>=0;
                 w = dbx1.frequency(getf);
                 w2 = dbx2.frequency(getf2);
                 switch rotation
                     case 'ww'
                         
                         %[W1,W2] = ndgrid(1:sum(getf),1:sum(getf2));      
                         [W1,W2] = ndgrid(1:sum(getf),1:bw2x:sum(getf2));      

                         W3 = round((w(W1)+w2(W2))./fstep3)+1;
                         
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
                 bspect = nan(size(W1));
                 I1=W1;I2=W2;I3=W3;
               
                      
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
%                   elseif  strcmpi(type,'bnb')
%                           cblrep = conj(blrep);
%                           blrep2 = dbx2.blrep(keept2,getf2,:);
%                   elseif  strcmpi(type,'bbb') || strcmpi(type,'nnn')
%                           cblrep = conj(blrep);
%                           blrep2 = blrep;
%                   end

                 bspect(inds) = sum(blrep(:,I1(inds),:).*blrep2(:,I2(inds),:).* cblrep(:,I3(inds),:));
                    
                 NRM = nan(size(I1));
                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds),:).*blrep2(:,I2(inds),:)).^2).* sum(abs(cblrep(:,I3(inds),:)).^2));
%                 NRM(inds) = sqrt(sum(abs(blrep(:,I1(inds),:)).^2).* sum(abs(cblrep(:,I3(inds),:).*blrep2(:,I2(inds),:)).^2));
             
end
   BICOH = bspect./NRM;
                     
   
   BICOH(isnan(BICOH))=0;
   
   
   