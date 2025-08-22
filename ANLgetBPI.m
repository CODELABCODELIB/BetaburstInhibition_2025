% ANLgetBBP: time locked beta burst probability
addpath(genpath('/home/wanw1/toolbox/eeglab2024'));

datapath = '/data1/projects/pi-ghosha2/Wenyu/study03/step03_burstdetection';
outpath  = '/data1/projects/pi-ghosha2/Wenyu/study03/step04_BBPI';
listing     = dir(datapath);
folderNames = {listing([listing.isdir]).name};
subID_list  = folderNames(~cellfun(@isempty,regexp(folderNames,'aaaa*')));
for s = 1: length(subID_list)
subID       = subID_list{1,s}

cd ([datapath '/' subID]);
setlist  = dir('1*.mat');
setname  = {setlist.name};

bbpi_all = [];

for i   = 1:length(setname)    
load(setname{1,i},'beta_burst_tap');
bbpi_all = cat(2,bbpi_all,beta_burst_tap);
end

bbpi_avg    = squeeze(nanmean(bbpi_all(:,:,:),2))'-median(squeeze(nanmean(bbpi_all(:,:,:),2))');

mkdir([outpath,'/' subID]);
cd([outpath,'/' subID]);

% plotting time-locked bbpi
load('/home/wanw1/toolbox/beta_burst/Orignalchanlocs.mat','Orignalchanlocs');
timepoints = [2000:200:4200]
figure('Units','normalized','Position',[0 0 1 1]);
for t = 1:length(timepoints)
subplot(2,6,t);
topoplot(bbpi_avg(timepoints(t),:), Orignalchanlocs,'maplimits',[1.1*min(bbpi_avg,[],'all') 1.1*max(bbpi_avg,[],'all')],'numcontour',0);
title([num2str(timepoints(t)-3000),'ms'])
end
%colorbar;
saveas(gca,'bbpi','png');

[min_peak,min_time]= min(bbpi_avg(2000:4000,:));
min_time           = min_time-1000;
% plotting negative peak
figure;topoplot(min_peak,Orignalchanlocs, 'numcontour',0);
colorbar;
saveas(gca,'bbpi_minpeak','png');
% plotting negpeak time
figure;topoplot(min_time,Orignalchanlocs, 'numcontour',0);
colorbar;
saveas(gca,'bbpi_minpeak_time','png');

[max_peak,max_time]= max(bbpi_avg(2000:4000,:));
max_time           = max_time-1000;
% plotting negative peak
figure;topoplot(max_peak,Orignalchanlocs, 'numcontour',0);
colorbar;
saveas(gca,'bbpi_maxpeak','png');
% plotting negpeak time
figure;topoplot(max_time,Orignalchanlocs, 'numcontour',0);
colorbar;
saveas(gca,'bbpi_maxpeak_time','png');

save('BBPI.mat','bbpi_avg','max_time','max_peak','min_time','min_peak');

close all
end