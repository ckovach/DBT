function [CH,tt,T] = choptf(trg,times,db,trref)

% [CH,tt] = choptf(trange,evnt_times,db,trref)
% Chop dbt data into epochs and normalize 
%
% trange - time epoch relative to...
% evnt_times - event markers
% db - dbt object
% trref - reference range by which to normalize. Normalization is currently
%         geometric mean amplitude envelope within the reference period applied to 
%         each corresponding trial epoch.
%         For no normalization leave empty (default)
% see also DBT

%C. Kovach 2014

if nargin < 4
    trref = [];
    nrm = 1;
end

adjust_phase = db.remodphase; %Correct phase for rounding error
                              %if phase information is preserved.
                              
[T,tt,Err]= chopper(trg,times,db.sampling_rate);

Err(T<1|T>length(db.time)) = 0;

T(T<1|T>length(db.time)) = length(db.time)+1;


if ~isempty(trref)

    Tref= chopper(trref,times,db.sampling_rate);
    Tref(Tref<1|Tref>length(db.time)) = length(db.time)+1;

else
    Tref = [];
end

if any(T(:)<1|T(:)>length(db.time)) || any(Tref(:)<1|Tref(:)>length(db.time))
    warning('%i time points are out of range of the input vector.\nCorresponding samples will be replaced with ones.',sum(T(:)<1|T(:)>length(db.time))+sum(Tref(:)<1|Tref(:)>length(db.time)))
end

x = ones(length(db.time)+1,1);
nch = size(db.blrep,3);
CH = zeros([length(tt),length(db.frequency),size(T,2), nch]);
for i = 1:nch
    for k = 1:length(db.frequency)

        x(1:length(db.time)) = db.blrep(:,k,i);

        if ~isempty(trref)
            %%% normalizing by mean envelope
%            nrm = repmat(mean(abs(x(Tref))),length(tt),1);
            nrm = repmat(exp(mean(log(abs(x(Tref))))),length(tt),1);
        end
        
        CH(:,k,:,i) = x(T)./nrm;
    end

end


if adjust_phase
    % This adjusts phase to compensate for rounding error
    adjuster  = exp(-1i*2*pi*repmat(permute(Err,[1 3 2]),[1,length(db.frequency),1]).*repmat(db.frequency - (~db.centerDC)*db.bandwidth*(1+db.shoulder)/2,[length(tt) 1 size(T,2)]));
    CH = CH.*repmat(adjuster,[1 1 1 size(CH,4)]);
end

