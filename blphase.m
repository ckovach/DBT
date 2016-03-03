function [dbout,PH,dbx]=blphase(x,fs,bw,varargin)

%
% [dbout,PH,dbx] = blphase(x,fs,bw)
%
% Band limited phase using iterative filtering 
% After each normalization the Hilbert transform is reapplied.
%
% dbout is a dbt object containing the normalized signal. H contains
% the reconstructed and amplitude normalized signal for each band at the
% original sampling rate.
%%
tol=1e-5;
shoulder = 1;
if isa(x,'dbt')
    fs = x.FSorig;
    bw = x.bandwidth;
    shoulder = x.shoulder;
    x = x.signal;
end
% DBT with hilbert equivalent trsfm
dbx = dbt(x,fs,bw,'centerdc',false,'shoulder',shoulder,varargin{:});

S = dbx.blrep./abs(dbx.blrep);

niter=0;
del=Inf;
dels=[];
while del > tol
    
     niter=niter+1;
    del
    
    F = fft(S);
    
    F(end/(2+2*dbx.fftpad)+1:end,:,:)=0;
    
    Snew=ifft(F);
    
    Snew=Snew./abs(Snew);
    del = sum(abs(Snew(:)-S(:)).^2)./sum(abs(S(:)).^2);
    
    S = Snew;
    
   % dels(end+1)=del;
end

dbout=dbx;
dbout.blrep=S;


if nargout > 1
    PH = zeros(length(x),length(dbout.frequency));
    for k = 1:length(dbout.frequency), 
        PH(:,k)=dbout.signal(k,1);
    end    
end
PH=PH./abs(PH);