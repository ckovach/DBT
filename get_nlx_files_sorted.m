function [fnames,fsizes] = get_nlx_files_sorted(get_path)


%%% Return nlx and nev files sorted numerically 

if nargin < 1
    get_path = false;
end



persistent path 

if ~get_path
    [fnames,pth] = uigetfile({'*.ncs;*.nev','NLX files (*.ncs,*.nev)';'*.nev','Event files only'},'Select NLX data file',path,'multiselect','on');
else
    pth = uigetdir(path,'Select NLX data file directory');
    fnames = cat(1,dir(fullfile(pth,'*.ncs')),dir(fullfile(pth,'*.nev')));
    fsizes=[fnames.bytes];
    fnames = {fnames.name};
end

if isnumeric(fnames)
    return
else
    path = pth;
end

if ~iscell(fnames)
    fnames = {fnames};
end

% sort numerically
decell= @(x)[x{:}];
reg = regexp(fnames,'(\d+)_*(\d*)','tokens','once');
reg(cellfun(@isempty,reg)) = {{'',''}};
chn = fliplr(arrayfun(@(x)str2double(['0',decell(x)]),cat(1,reg{:})));
[~,srti] = sortrows(chn);

fnames = fullfile(path,fnames(srti));

if get_path
    fsizes = fsizes(srti);
else
    fsizes = nan;
end