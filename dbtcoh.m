function [coh,csp,w,dbs,trf] =dbtcoh(x,y,varargin)


% Simple function to compute coherence values with dbt transforms.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


nx = size(x,2);    

if isscalar(y)
    varargin = [{y},varargin];
    y=[];
    ny = nx;
elseif ~isempty(y) 
    dby = dbt(y,varargin{:});
    ny = size(y,2);
else
    ny = nx;
end

dbx = dbt(x,varargin{:});

w = dbx.frequency;
csp = zeros(nx,ny,length(w));
coh = csp;
if nargout > 4
    trf = csp;
end
for i = 1:length(dbx.frequency)
        
        blx = squeeze(dbx.blrep(:,i,:));
        if isempty(y)
            bly = blx;
        else
            bly = squeeze(dby.blrep(:,i,:));
        end
        csp(:,:,i) = blx'*bly;
        if isempty(y)
             coh(:,:,i) = diag(diag(csp(:,:,i).^-.5))*csp(:,:,i)*diag(diag(csp(:,:,i)).^-.5);
        else          
           coh(:,:,i) = diag(sum(abs(blx).^2).^-.5)*csp(:,:,i)*diag(sum(abs(bly).^2).^-.5);
        end
        
    if nargout > 4
           trf(:,:,i) = diag(sum(abs(blx).^2))\(blx'*bly);
          
    end
end
    
if nargout > 3
    dbs = [dbx,dby];
end



