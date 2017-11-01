
function [DSS,SSR,DTF,dbxpast,dbxfut]  = dbt_causality(dbx,varargin)

% This is an experimental application of causal modeling using a crude
% but efficient factorization.



make_plot = true;

if ~dbx.remodphase
    error('dbx.remodphase must be TRUE');
end

B = dbx.blrep;

FB = ifft(B,[],2);

nfr = length(dbx.frequency);
w = ifftshift((0:nfr-1)-floor(.5*nfr))/nfr;

Bpast = FB;
Bpast(:,w>=0,:) = 0;

Bfut = FB;
Bfut(:,w<=0,:) = 0;

dbxpast = dbx;
dbxpast.blrep = fft(Bpast,[],2);
dbxfut = dbx;
dbxfut.blrep = fft(Bfut,[],2);

DTF = dbtcoh(dbxpast,dbxfut);


for k = 1:length(dbx.frequency);
    
    Y = squeeze(dbxfut.blrep(:,k,:));
    X = squeeze(dbxpast.blrep(:,k,:));
    ssY = sum(abs(Y).^2);
    ssX = sum(abs(X).^2);
    
    npar = size(X,2);
    for kk = 1:npar+1;
        getw = (1:npar)~=kk-1;
        [glm(k,kk).b,glm(k,kk).dev,glm(k,kk).pval,glm(k,kk).iXX,glm(k,kk).sigma,glm(k,kk).res,glm(k,kk).Yfit] = complexglm(Y,X(:,getw),'intercept',true);

        ssR = diag(ssX)*abs(glm(k,kk).b(1:size(ssX,2),:)).^2*diag(ssY.^-1);
        glm(k,kk).ssR =ssR;
        glm(k,kk).ssRes = sum(abs(glm(k,kk).res).^2);
        
        glm(k,kk).dSSR = (glm(k,kk).ssRes-glm(k,1).ssRes)./ssY;
    end%     if any(ssR(:)>1)
%         keyboard
%     end
    k 
end
for k =1:npar
    DSS(:,:,k) = cat(1,glm(:,k+1).dSSR);
    DSS(:,k,k) = 0;
end
SSR = cat(3,glm(:,1).ssR);

if make_plot
    figure
    pln = 1;
    yl = minmax(DSS(:));
    for k = 1:npar
        for kk = 1:npar
            if k~=kk
                
                subplot(npar,npar,pln)
                plot(dbx.frequency,squeeze(DSS(:,kk,k)))
                title(sprintf('%i -> %i',k,kk))
                ylim(yl)
            end
            pln = pln+1;
        end
    end
end
