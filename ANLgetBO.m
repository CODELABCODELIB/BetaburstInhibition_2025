% ANLgetBO: burst occupancy
addpath(genpath('/home/wanw1/toolbox/eeglab2024'));
addpath(genpath('/home/wanw1/toolbox/morans_i-master'));

datapath = '/data1/projects/pi-ghosha2/Wenyu/study03/step03_burstdetection';
outpath  = '/data1/projects/pi-ghosha2/Wenyu/study03/step04_BBP';
listing     = dir(datapath);
folderNames = {listing([listing.isdir]).name};
subID_list  = folderNames(~cellfun(@isempty,regexp(folderNames,'aaaa*')));
for s = 1: length(subID_list)
subID       = subID_list{1,s}
cd ([datapath '/' subID]);
setlist  = dir('1*.mat');
setname  = {setlist.name};

beta_burst_dur_all = zeros(1,62);
recording_dur_all  = zeros(1,62);
for i   = 1:length(setname)    
load(setname{1,i},'beta_burst_dur','recording_dur');
beta_burst_dur_all = beta_burst_dur_all+beta_burst_dur;
recording_dur_all  = recording_dur_all+recording_dur;
end

bbp    = beta_burst_dur_all./recording_dur_all;

mkdir([outpath,'/' subID]);
cd([outpath,'/' subID]);
save('BBP.mat','bbp');

% plotting
load('/home/wanw1/toolbox/beta_burst/Orignalchanlocs.mat','Orignalchanlocs');
figure;topoplot(bbp, Orignalchanlocs,'maplimits',[0.9*min(bbp) 1.1*max(bbp)],'numcontour',0);
colorbar;
saveas(gca,'bbp','png');

end