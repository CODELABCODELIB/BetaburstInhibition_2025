% ANLgetBBP
addpath(genpath('/home/wanw1/toolbox/eeglab2024'));
addpath(genpath('/home/wanw1/toolbox/morans_i-master'));

datapath = '/data1/projects/pi-ghosha2/Wenyu/study03/step03_burstdetection';
outpath  = '/data1/projects/pi-ghosha2/Wenyu/study03/step04_BBR';
listing     = dir(datapath);
folderNames = {listing([listing.isdir]).name};
subID_list  = folderNames(~cellfun(@isempty,regexp(folderNames,'aaaa*')));
for s = 1: length(subID_list)
subID       = subID_list{1,s}
%subID    = 'aaaaaaeF';
cd ([datapath '/' subID]);
setlist  = dir('1*.mat');
setname  = {setlist.name};

bb_tapnum_all       = zeros(1,62);
bb_dur_all          = zeros(1,62);
rec_dur_all         = zeros(1,62);
tapnum_all          = zeros(1,62);
for i   = 1:length(setname)    
load(setname{1,i},'beta_burst_dur','tapnum','recording_dur','iti');
bb_tapnum_all         = bb_tapnum_all+tapnum;
bb_dur_all            = bb_dur_all+beta_burst_dur;
rec_dur_all           = rec_dur_all+recording_dur;
tapnum_all            = repmat(size(iti,2),1,62)+tapnum_all;
end

bbr                = tapnum_all./bb_dur_all;
nbr                = (tapnum_all-bb_tapnum_all)./(rec_dur_all-bb_dur_all);

gap_br             = nbr-bbr;

mkdir([outpath,'/' subID]);
cd([outpath,'/' subID]);
save('BBR.mat','bbr','nbr','gap_br');

% plotting
load('/home/wanw1/toolbox/beta_burst/Orignalchanlocs.mat','Orignalchanlocs');
figure;topoplot(bbr, Orignalchanlocs,'maplimits',[min(bbr) max(bbr)],'numcontour',0);
colorbar;
saveas(gca,'bbr','png');

figure;topoplot(nbr, Orignalchanlocs,'maplimits',[min(nbr) max(nbr)],'numcontour',0);
colorbar;
saveas(gca,'nbr','png');

figure;topoplot(nbr-bbr, Orignalchanlocs,'numcontour',0);%'maplimits',[0 1.1*max(nbr-bbr)],
colorbar;
saveas(gca,'nbr-bbr','png');


end