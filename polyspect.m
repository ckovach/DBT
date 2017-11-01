function [PCOH,w,pspect,NRM,dbx1] = polyspect(x,fs,bw,ord,varargin)

% Computes the nth order polyspectrum.
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

maxfreq=Inf;
warning('This script is under development and is not stable. Do not use it unless you know what you''re doing')
symmetric= false;
if nargin < 3 || isempty(bw)
    bw = 1;      
end

if nargin < 4 || isempty(ord)
   ord = 3; 
end

upsampfx = 1;
                
                
getfreq = min(maxfreq,fs/2);
upsamptx = 1;
dbx1 = dbt(x,fs,bw,'upsampleFx',upsampfx,'remodphase',false,'upsampleTx',upsamptx-1,'lowpass',getfreq);
getf = dbx1.frequency<=getfreq & dbx1.frequency>=0;
Ws = repmat({[]},1,ord-1);
%      switch rotation
%          case 'ww'

[Ws{:}] = ndgrid(1:sum(getf));      
Ws{ord} = sum(cat(ord,Ws{:}),ord)-ord+2;
             

%          case 'vv'
%              [V1,V2] = ndgrid(1:sum(getf));        
%             W1=(V1-V2)+1;
%             W2=(V1+V2)-1;
%             W3 = 2*V1-1;
% 
%      end
     %  inds = W3<=length(dbx.frequency);
 inds = true;
 for k = 1:ord-1
     inds= inds & (Ws{k}>0 & Ws{k} <=sum(getf));
     if k<ord-1 && ~symmetric
         inds = inds& Ws{k}>=Ws{k+1};
     end
 end
 inds= inds & Ws{ord}<=sum(dbx1.frequency>=0);


 w = dbx1.frequency(getf);
 blrep = dbx1.blrep(:,dbx1.frequency>=0,:);
 pspect = nan([size(Ws{1}),size(blrep,3)]);

 B = ones(size(blrep,1),sum(inds(:)),size(blrep,3));
 for k = 1:ord-1
     B = B.*blrep(:,Ws{k}(inds),:);
 end
 C = conj(blrep(:,Ws{ord}(inds),:));
 pspect(inds) = sum(B.*C);
        
 NRM = nan(size(Ws{1}));
 NRM(inds) = sqrt(sum(abs(B).^2).*sum(abs(C).^2));


PCOH = pspect./NRM;


PCOH(isnan(PCOH))=0;

sum(inds(:))

