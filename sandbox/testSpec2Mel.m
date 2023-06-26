addpath ../utils/
addpath ../tools/
addpath ../scripts/

par=projParam();

spfile=[par.pth.speechpth '/timit/train/dr8/mbsb0/si723.wav'];
[s,fs]=v_readsph(spfile,'wt');
tt=(0:length(s)-1)/fs;

%% 1) Get frames
framesEpoch=getFrames(s,fs,par.db,2);
framesFixed=getFrames(s,fs,par.db,1);

%% 2) Calculate STFTs
[stftFixed,faxFixed]=gs_stft(s,framesFixed,par.db.nfft);
[stftEpoch,faxEpoch]=gs_stft(s,framesEpoch,par.db.nfft);
[stftEpoch_nzp,framelen]=gs_stft_nzp(s,framesEpoch,par.db.nfft);

%% 3) Convert spectrograms to Mel energy - grams
numMels=22;
preserveDC=0;
stftmFixed=spec2mel(stftFixed,fs,numMels,preserveDC);
stftmEpoch=spec2mel(stftEpoch,fs,numMels,preserveDC);
stftm_nzp=spec2mel_nzp(stftEpoch_nzp,fs,numMels,framelen);

%% 4) Convert Mel Energy-grams to spectrograms
stftrFixed=mel2spec(stftmFixed,fs,par.db.nfft,preserveDC,angle(stftFixed));
stftrEpoch=mel2spec(stftmEpoch,fs,par.db.nfft,preserveDC,angle(stftEpoch));
stftrEpoch_nzp=mel2spec_nzp(stftm_nzp,fs,par.db.nfft,framelen,angle(stftEpoch_nzp));

%% 5) Invert spectrograms to time domain signal
vFixed=gs_istft(stftrFixed,framesFixed);
vEpoch=gs_istft(stftrEpoch,framesEpoch);
vEpoch_nzp=gs_istft_nzp(stftrEpoch_nzp,framesEpoch);

%% 6) Calculate SNR
SNRfixed=calculateSNR(s,vFixed);
SNRepoch=calculateSNR(s,vEpoch);
SNRepoch_nzp=calculateSNR(s,vEpoch_nzp);

disp([SNRfixed, SNRepoch, SNRepoch_nzp]);


