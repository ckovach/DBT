function [CH,tt] = choptf(trg,times,db,trref)

% [CH,tt] = choptf(trange,evnt_times,db,trref)
% Chop dbt data into epochs and normalize 
%
% trange - time epoch relative to...
% evnt_times - event markers
% db - dbt object
% trref - reference range by which to normalize. Normalization is currently
%         mean amplitude envelope within the reference period applied to 
%         each corresponding trial epoch.
%         For no normalization leave empty (default)
% see also DBT

%C. Kovach 2014

if nargin < 4
    trref = [];
end

[T,tt]= chopper(trg,times,db.sampling_rate);
T(T<1|T>length(db.time)) = length(db.time)+1;

if ~isempty(trref)

    Tref= chopper(trref,times,db.sampling_rate);
    Tref(Tref<1|Tref>length(db.time)) = length(db.time)+1;

end
x = ones(length(db.time)+1,1);
CH = zeros([length(tt),length(db.frequency),size(T,2)]);
for k = 1:length(db.frequency)
    
    x(1:length(db.time)) = db.blrep(:,k);
    
    if ~isempty(trref)
        %%% normalizing by mean envelope
        nrm = repmat(mean(abs(x(Tref))),length(tt),1);
    end
    CH(:,k,:) = x(T)./nrm;
            

end
