function [coh,csp,w,dbs,trf] =dbtcoh(x,y,varargin)


% Simple function to compute coherence values with dbt transforms.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

dbx = dbt(x,varargin{:});

dby = dbt(y,varargin{:});

w = dbx.frequency;
csp = zeros(size(x,2),size(y,2),length(w));
coh = csp;
if nargout > 4
    trf = csp;
end
for i = 1:length(dbx.frequency)
        
        blx = squeeze(dbx.blrep(:,i,:));
        bly = squeeze(dby.blrep(:,i,:));
        
        csp(:,:,i) = blx'*bly;
        coh(:,:,i) = diag(diag(csp(:,:,i).^-.5))*csp(:,:,i)*diag(diag(csp(:,:,i)).^-.5);
    if nargout > 4
           trf(:,:,i) = diag(sum(abs(blx).^2))\(blx'*bly);
          
    end
end
    
if nargout > 3
    dbs = [dbx,dby];
end



