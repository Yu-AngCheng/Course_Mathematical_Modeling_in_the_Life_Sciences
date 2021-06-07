for k = 1:39
    dat=alldat{k};
    brain_groups = {{'VISa', 'VISam', 'VISl', 'VISp', 'VISpm', 'VISrl'},
        {'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'POL', 'PT', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM'}
        {'CA', 'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST'},
        {'ACA', 'AUD', 'COA', 'DP', 'ILA', 'MOp', 'MOs', 'OLF', 'ORB', 'ORBm', 'PIR', 'PL', 'SSp', 'SSs', 'RSP',' TT'}
        {'APN', 'IC', 'MB', 'MRN', 'NB', 'PAG', 'RN', 'SCs', 'SCm', 'SCig', 'SCsg', 'ZI'},
        {'ACB', 'CP', 'GPe', 'LS', 'LSc', 'LSr', 'MS', 'OT', 'SNr', 'SI'},
        {'BLA', 'BMA', 'EP', 'EPd', 'MEA'}
        };
    t=size(dat.spks);
    nArea=t(1);
    nTrial=t(2);
    nBin=t(3);
    group=ones(nArea,1)*7;
    for i=1:nArea
        name=dat.brain_area(i,:);
        name(isspace(name))=[];
        for j=1:6
            if ismember(name,brain_groups{j})
                group(i)=j;
            end
        end
    end
    groupSpk=zeros(7,nTrial,nBin);
    t=size(dat.spks_passive);
    groupSpkPassive=zeros(7,t(2),nBin);
    for i=1:7
        groupSpk(i,:,:)=mean(dat.spks(group==i,:,:));
        groupSpkPassive(i,:,:)=mean(dat.spks_passive(group==i,:,:));
    end
    dat.group=group;
    dat.groupSpk=groupSpk;
    dat.groupSpkPassive=groupSpkPassive;
    alldat{k}=dat;
end
         