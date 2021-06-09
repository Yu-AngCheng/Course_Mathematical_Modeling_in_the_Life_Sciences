%% load data
if exist('Steinmetz_main_braingroup.mat')
    error('Already preprocessed')
end
clear
load Steinmetz_main.mat
%% set the parameters
movingWindow = 5;
movingmethod = 'movmean';
brain_groups = {
    'VISa', 'VISam', 'VISl', 'VISp', 'VISpm', 'VISrl',...
    'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'POL', 'PT', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM',...
    'CA', 'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST',...
    'ACA', 'AUD', 'COA', 'DP', 'ILA', 'MOp', 'MOs', 'OLF', 'ORB', 'ORBm', 'PIR', 'PL', 'SSp', 'SSs', 'RSP',' TT',...
    'APN', 'IC', 'MB', 'MRN', 'NB', 'PAG', 'RN', 'SCs', 'SCm', 'SCig', 'SCsg', 'ZI',...
    'ACB', 'CP', 'GPe', 'LS', 'LSc', 'LSr', 'MS', 'OT', 'SNr', 'SI',...
    'BLA', 'BMA', 'EP', 'EPd', 'MEA'
    };
nGroup = length(brain_groups);
%% calculate the firing rate by every brain group
for k = 1:39

    dat = alldat{k};
    
    dat.spks = permute(dat.spks,[2,3,1]);    
    tmp = size(dat.spks);
    nArea = tmp(3);nTrial = tmp(1);nBin = tmp(2);
    
    group = ones(nArea,1)*(nGroup+1);
    for i = 1:nArea
        name = dat.brain_area(i,:);name(isspace(name)) = [];
        [isin,idx] = ismember(name,brain_groups);
        if isin ~= 0
            group(i) = idx;
        end
    end
    
    groupSpk = NaN(nTrial,nBin,nGroup+1);
    for i = 1:(nGroup+1)
        groupSpk(:,:,i) = nanmean(dat.spks(:,:,group==i),3);
    end
    
    dat.group = group;
    dat.groupSpk = groupSpk;
    dat.groupFiringRate = smoothdata(dat.groupSpk,2,movingmethod,movingWindow);
    dat.session = repmat(k,nTrial,1);
    
    alldat{k} = dat;
    clearvars -except alldat movingWindow brain_groups nGroup k movingmethod
end
save Steinmetz_main_braingroup.mat alldat -v7.3