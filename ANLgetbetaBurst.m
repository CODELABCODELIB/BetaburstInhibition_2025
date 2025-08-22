function ANLgetbetaBurst(inpath, outpath,subID)
% Wenyu wan, Leiden univerisity, 09/2024, edited: 15/05/2025
cd ([inpath subID]);
setlist      = dir('1*.set');
setname      = {setlist.name};

% inter-beta burst onset interval JID
for s = 1:length(setname)
outdir         = strcat(outpath,subID);
mkdir(outdir)

load('/home/wanw1/toolbox/beta_burst/Orignalchanlocs.mat');
% load dataset/d
cd ([inpath subID]);
EEG           = pop_loadset(strcat(inpath, subID,'/', setname{1,s}));

% remove artifact IC
outdir_tmp      = [pwd filesep 'amicaout'];
cd(outdir_tmp);
modout          = loadmodout15(outdir_tmp);
model_index     = 1;
EEG.icawinv     = modout.A(:,:,model_index);
EEG.icaweights  = modout.W(:,:,model_index);
EEG.icasphere   = modout.S;
EEG             = eeg_checkset(EEG);
EEG             = pop_iclabel(EEG,'default');
EEG             = pop_icflag(EEG,[NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
EEG             = pop_subcomp(EEG,[],0,0);

% remove eye channels
EEG             = pop_select(EEG, 'nochannel',{'E5' 'E64'});

% beta-band filtering
EEG             = pop_eegfiltnew(EEG, [],45);
EEG             = pop_eegfiltnew(EEG, 1,[]);

% remove unrelevant channels
%EEG            = pop_rejchan(EEG, 'threshold', 3, 'norm', 'on', 'measure', 'spec','freqrange',[13 30]);
EEG             = clean_channels(EEG,0.85);

% interpolate
EEG             = pop_interp(EEG,Orignalchanlocs,'spherical');

% reject inactive data segment
EEG             = rejectsegment(EEG,0,0);

% epoch
INEEG           = pop_epoch(EEG,{'Phone'},[-3 3],'epochinfo','no');
INEEG           = selectepoch(INEEG);

beta_burst_tap       = nan(62,length(INEEG.modepoch),6001);
iti                  = nan(62,length(INEEG.modepoch));
burst_rate_iti       = nan(62,length(INEEG.modepoch));

% Define parameters
for c= 1:62       % Channel to analyze
fre_band          = [1:1:45];            % Beta-band frequency range (Hz)
beta_band         = [13:1:30];
fs                = 1000;
m                 = 7; 
time              = -1:1/fs:1;
%Create time-frequency power matrix
wavelet_conv_data = zeros(length(fre_band),EEG.pnts);
% from paper
for f_idx = 1:length(fre_band)
    % freqs and cycles
    freq = fre_band(f_idx);
    sigma = m/(2*pi*freq);
    % sine and gaussian
    sine_wave = exp(1i*2*pi*freq.*time);
    gaussian_win = exp(-time.^2./(2*sigma^2));
    % normalization factor
    normalization_factor = 1 / (sigma * sqrt(2* pi));
    % make wavelet
    wavelet = normalization_factor .* sine_wave .* gaussian_win;
    halfwaveletsize = ceil(length(wavelet)/2); % half of the wavelet size
    % convolve with data
    n_conv = length(wavelet) + EEG.pnts - 1; % compute Gaussian
    % fft
    fft_w = fft(wavelet,n_conv);
    fft_e = fft(EEG.data(c,:),n_conv);
    ift   = ifft(fft_e.*fft_w,n_conv);
    wavelet_conv_data(f_idx,:) = abs(ift(halfwaveletsize:end-halfwaveletsize+1)).^2;
end

detected_tf      = zeros(size(wavelet_conv_data,1),size(wavelet_conv_data,2));
% power threshold
burst_cutoff      =  6*nanmedian(wavelet_conv_data');
% detect beta-burst based on power threshold
for f_idx = 1:length(beta_band)
    detected_tf(f_idx,:) = wavelet_conv_data(f_idx,:)>burst_cutoff(f_idx);
end
% futher check based on duration burst

detected_tf_new         = zeros(size(detected_tf,1),size(detected_tf,2));
for f_idx = 1:length(fre_band)
    detected_tf(f_idx,:) = wavelet_conv_data(f_idx,:)>burst_cutoff(f_idx);
    
    betaBurstInds = SplitVec(find(detected_tf(f_idx,:)),'consecutive');
    segL          = cellfun('length',betaBurstInds); % Find burst lengths
    burstSelInds  = segL>((1000/fre_band(f_idx))*2); % Select bursts with above min length
    burstSelIndsout = betaBurstInds(burstSelInds);
     for i         = 1:length(burstSelIndsout)
         detected_tf_new(f_idx,burstSelIndsout{1,i}) = 1;   
     end
end

detected_idx           = (nansum(detected_tf_new(beta_band,:),1)>0);
freq_prob(c,:)         = nanmean(detected_tf_new(beta_band,:),2);
% [M,I]                  = max(freq_prob(c,:));
% freq_M(c)              = M;
% freq_I(c)              = beta_band(I);

for e     = 2: length(INEEG.modepoch)
     try
         beta_burst_tap(c,e,:)    = detected_idx(INEEG.urevent(INEEG.modepoch(e).eventurevent).latency-3000:INEEG.urevent(INEEG.modepoch(e).eventurevent).latency+3000);
         %tf_tap(c,e,:,:)          = wavelet_conv_data(beta_band,INEEG.urevent(INEEG.modepoch(e).eventurevent).latency-3000:INEEG.urevent(INEEG.modepoch(e).eventurevent).latency+3000);
         iti(c,e)               = (INEEG.urevent(INEEG.modepoch(e).eventurevent).latency-INEEG.urevent(INEEG.modepoch(e-1).eventurevent).latency);
         burst_rate_iti(c,e)      = nansum(detected_idx(INEEG.urevent(INEEG.modepoch(e-1).eventurevent).latency:INEEG.urevent(INEEG.modepoch(e).eventurevent).latency))/iti(c,e);
         tap_idx(e)               = INEEG.urevent(INEEG.modepoch(e).eventurevent).latency;
         %iti_1(e)                = (INEEG.urevent(INEEG.modepoch(e+1).eventurevent).latency-INEEG.urevent(INEEG.modepoch(e).eventurevent).latency);
     end
 end

 
 % burst rate: burst duration/whole duration
 % burst_rate_whole(c)   = sum(detected_idx>0)/length(detected_idx);
 beta_burst_dur(c)            = nansum(detected_idx>0);
 beta_burst_dur_allfre(c,:)   = nansum(detected_tf_new(beta_band,:)>0,2);
 recording_dur(c)             = length(detected_idx);
 
 % calculate tap num during beta burst
 betaburst_idx           = find(detected_idx>0);
 tapnum(c)               = nansum(ismember(tap_idx,betaburst_idx));

end
cd (outdir);
filename                 = erase(setname{1,s},'.set');
save(strcat(filename,'.mat'),'burst_rate_iti','iti','beta_burst_dur','beta_burst_dur_allfre','recording_dur','beta_burst_tap','freq_prob','tapnum','-v7.3');
end

end
