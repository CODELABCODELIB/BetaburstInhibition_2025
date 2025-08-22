% ANLgetBBTI: burst behavior timing index
addpath(genpath('/home/wanw1/toolbox/eeglab2024'));
% addpath(genpath('/home/wanw1/toolbox/morans_i-master'));

datapath = '/data1/projects/pi-ghosha2/Wenyu/study03/step03_burstdetection';
outpath  = '/data1/projects/pi-ghosha2/Wenyu/study03/step05_ITIlog';
listing     = dir(datapath);
folderNames = {listing([listing.isdir]).name};
subID_list  = folderNames(~cellfun(@isempty,regexp(folderNames,'aaaa*')));
for s = 1: length(subID_list)
subID       = subID_list{1,s}
cd ([datapath '/' subID]);
setlist  = dir('1*.mat');
setname  = {setlist.name};

bbiti_all = [];
iti_all   = [];
for i     = 1:length(setname)    
load(setname{1,i},'burst_rate_iti','iti');
bbiti_all  = cat(2,bbiti_all,burst_rate_iti);
iti_all    = cat(2,iti_all,iti);
end


mkdir([outpath,'/' subID]);
cd([outpath,'/' subID]);

for c = 1:62
idx_burst   = bbiti_all(c,:)>0;
iti_burst   = iti_all(1,idx_burst);
iti_noburst = iti_all(1,~idx_burst);
[h(c),p(c),~,tstats] = ttest2(log10(iti_burst),log10(iti_noburst));
stats(c)    = tstats.tstat;
iti_gap(c)  = nanmean(log10(iti_burst))-nanmean(log10(iti_noburst));
iti_burstavg(c)   = nanmean(log10(iti_burst));
iti_noburstavg(c) = nanmean(log10(iti_noburst));
end


load('/home/wanw1/toolbox/beta_burst/Orignalchanlocs.mat','Orignalchanlocs');
% plotting difference iti distribution between two groups
figure;topoplot(iti_gap, Orignalchanlocs,'maplimits',[min(iti_gap)*0.9 max(iti_gap)*1.1],'numcontour',0);
colorbar;
saveas(gca,'itigap_distribution','png');


% plotting statistic
figure;topoplot(stats, Orignalchanlocs,'maplimits',[min(stats)-0.5 max(stats)+0.5],'numcontour',0);
colorbar;
saveas(gca,'ttest_distribution','png');

save('BBITI.mat','p','h','stats','bbiti_all','iti_all','iti_gap','iti_burstavg','iti_noburstavg');

iti_allsub(s,1) = nanmedian(iti_all(1,:));
iti_allsub(s,2) = iqr(iti_all(1,:));
close all
end