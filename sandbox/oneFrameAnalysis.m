addpath ../utils/
addpath ../tools/
addpath ../scripts/

par=projParam();

spfile=[par.pth.speechpth '/timit/train/dr8/mbsb0/si723.wav'];
[s,fs]=v_readsph(spfile,'wt');
tt=(0:length(s)-1)/fs;
wl=round(25e-3*fs);  % 25 ms in samples

framesEpoch=getFrames(s,fs,par.db,2);
framesFixed=getFrames(s,fs,par.db,1);

[~,ixEpoch]=min(abs(framesEpoch/fs-1.51)); %Pick the frame at ca 1.51 s
[~,ixFixed]=min(abs(framesFixed/fs-1.51)); %Pick the frame at ca 1.51 s
nnsegEpoch=framesEpoch(ixEpoch):framesEpoch(ixEpoch+1)-1; %Sample indices of the chosen 
nnsegFixed=framesFixed(ixFixed):(framesFixed(ixFixed+1)-1);

figure(1); clf;
stem(framesEpoch/fs,max(s)*ones(size(framesEpoch))); hold on
plot(tt,s);
plot(tt(nnsegFixed),s(nnsegFixed),'k');%,'LineWidth',1);
plot(tt(nnsegEpoch),s(nnsegEpoch),'--','Color',0.75*[1 1 1],'LineWidth',1.2);
xlim([1.49 1.54]);
xlabel('Time [s]')
ylabel('Normalised speech signal amplitude');

%% Fourier analysis
%
% Use gs_stft and gs_stft_noPadding to calculate the spectrograms.  Pick
% out the frame of interest and plot power spectrum, phase and group-delay.

[stftEpoch,fax]=gs_stft(s,framesEpoch,par.db.nfft);
[stftFixed,fax]=gs_stft(s,framesFixed,par.db.nfft);
[stft,framelen]=gs_stft_noPadding(s,framesEpoch,par.db.nfft);
ixFlen=framelen(ixEpoch)/2;
faxNoPadd=(0:ixFlen)*fs/ixFlen*0.5;
nfftp=1+floor(ixFlen/2); 

figure(2); clf;
ax2(1)=subplot(311);
plot(fax*fs,20*log10(abs(stftFixed(1:length(fax),ixFixed)))); hold on;
plot(fax*fs,20*log10(abs(stftEpoch(1:length(fax),ixEpoch))));
plot(faxNoPadd,20*log10(abs(stft(1:length(faxNoPadd),ixEpoch))));
ylabel('Power [dB]')
ax2(2)=subplot(312); 
plot(fax*fs,unwrap(angle(stftFixed(1:length(fax),ixFixed)))); hold on;
plot(fax*fs,unwrap(angle(stftEpoch(1:length(fax),ixEpoch))));
plot(faxNoPadd,unwrap(angle(stft(1:length(faxNoPadd),ixEpoch))));
ylabel('Unwrapped angle [rad]')
ax2(3)=subplot(313); 
plot(fax*fs,[0; diff(unwrap(angle(stftFixed(1:length(fax),ixFixed))))]*nfftp); hold on;
plot(fax*fs,[0; diff(unwrap(angle(stftEpoch(1:length(fax),ixEpoch))))]*nfftp);
plot(faxNoPadd,[0; diff(unwrap(angle(stft(1:length(faxNoPadd),ixEpoch))))]*ixFlen);
ylabel('Group Delay')
xlabel('Frequency [Hz]')
linkaxes(ax2,'x')

%% Mel analysis 1
%
% Plot the melbank and the spectrum together.  

numMels=22;
[mbmPadd, fxmbmPadd]=v_filtbankm(numMels,par.db.nfft,fs,0,fs/2,'m');
[mbmNoPadd, fxmbmNoPadd]=v_filtbankm(numMels,framelen(ixEpoch),fs,0,fs/2,'m');
%fxmbmPadd and fxmbmNoPadd should be the same since numMels and fs are the
%same.

figure(3); clf;
ax3(1)=subplot(311);
plot(fax*fs,mbmPadd'); hold on;
plot(fax*fs,(abs(stftFixed(1:length(fax),ixFixed))),'k'); 
ylabel('Magnitude')
ax3(2)=subplot(312);
plot(fax*fs,mbmPadd'); hold on;
plot(fax*fs,(abs(stftEpoch(1:length(fax),ixEpoch))),'k');
ax3(3)=subplot(313);
plot(faxNoPadd',mbmNoPadd'); hold on;
plot(faxNoPadd',(abs(stft(1:length(faxNoPadd),ixEpoch))),'k');
linkaxes(ax3','x')

imbmPadd=pinv(full(mbmPadd));
imbmNoPadd=pinv(full(mbmNoPadd));

% Energy in the mel-filters
nrgmbm=sqrt(sum(full(mbmPadd')));
nrgmbmNoPadd=sqrt(sum(full(mbmNoPadd')));

%% Mel analysis 2

nfftp=length(fax);
[sdftrFixed, sdftmFixed]=melAndReconstruct(stftFixed,fs,par.db.nfft,nfftp,numMels,par.db.fbankmethod,par.db.preserveDC);
[sdftrEpoch, sdftmEpoch]=melAndReconstruct(stftEpoch,fs,par.db.nfft,nfftp,numMels,par.db.fbankmethod,par.db.preserveDC);
stftm=spec2melNoPadding(stft,fs,numMels,framelen);

% We compare the Mel-energies and notice how the nonpadded epoch mel
% cepstrum is similar to the padded one (with a ratio of nfftp/ixFlen)
figure(4); clf;
plot((0:numMels-1),sdftmFixed(:,ixFixed)); hold on;
plot((0:numMels-1),sdftmEpoch(:,ixEpoch));
plot((0:numMels-1),nfftp/ixFlen*stftm(:,ixEpoch));
ylabel('Mel filter energy')
xlabel('Mel index')

figure(41); clf;
plot(fxmbmPadd,sdftmFixed(:,ixFixed)); hold on;
plot(fxmbmPadd,sdftmEpoch(:,ixEpoch));
plot(fxmbmNoPadd,nfftp/ixFlen*stftm(:,ixEpoch));
xlabel('Frequency [Hz]')
ylabel('Mel filter energy')

% If we can plot the mel filter bank values like this, why not compare them
% to the spectrum:

figure(42); clf;
ax42(1)=subplot(311);
plot(fax*fs,(abs(stftFixed(1:length(fax),ixFixed))),'k'); hold on;
plot(fxmbmPadd,sdftmFixed(:,ixFixed)'./nrgmbm,'b*');
ylabel('Magnitude')
ax42(2)=subplot(312);
plot(fax*fs,(abs(stftEpoch(1:length(fax),ixEpoch))),'k'); hold on;
plot(fxmbmPadd,sdftmEpoch(:,ixEpoch)'./nrgmbm,'b*');
ax42(3)=subplot(313);
plot(faxNoPadd',nfftp/ixFlen*(abs(stft(1:length(faxNoPadd),ixEpoch))),'k'); hold on;
plot(fxmbmNoPadd,nfftp/ixFlen*stftm(:,ixEpoch)'./nrgmbmNoPadd,'b*');
linkaxes(ax42','x')


%% Reconstructed spectrum
%
% melAndReconstruct give us the reconstructed spectrum using pinv on mbm.
% mel2specNoPadding also implments pinv to get a reconstructed spectrum.
% In the script we also experiment with using interpolation.

% Running mel2specNoPadding so we get the pinv version of the reconstructed
% "nopadding" spectrum
stftangle=angle(stft);
stftr=mel2specNoPadding(stftm,fs,par.db.nfft,framelen,stftangle);

% Use interpolation directly just on that single frame.
%Vq = interp1(X,V,Xq);
imet='linear';
XfixedR = interp1(fxmbmPadd,(sdftmFixed(:,ixFixed)./nrgmbm),fax*fs,imet,'extrap');
XepochR = interp1(fxmbmPadd,(sdftmEpoch(:,ixFixed)./nrgmbm),fax*fs,imet,'extrap');
XenopadR = interp1(fxmbmNoPadd,stftm(:,ixFixed),faxNoPadd,imet,'extrap');

figure(6); clf;
ax5(1)=subplot(311);
plot(fax*fs,20*log10(abs(stftFixed(1:length(fax),ixFixed)))); hold on;
plot(fax*fs,20*log10(abs(sdftrFixed(1:length(fax),ixFixed))));
plot(fax*fs,20*log10(max(XfixedR(1:length(fax)),0)));
ylabel('Power [dB]')
ax5(2)=subplot(312);
plot(fax*fs,20*log10(abs(stftEpoch(1:length(fax),ixEpoch)))); hold on;
plot(fax*fs,20*log10(abs(sdftrEpoch(1:length(fax),ixEpoch))));
plot(fax*fs,20*log10((XepochR(1:length(fax)))));
ylabel('Power [dB]')
ax5(3)=subplot(313);
plot(faxNoPadd,20*log10(abs(stft(1:length(faxNoPadd),ixEpoch)))); hold on;
plot(faxNoPadd,20*log10(abs(stftr(1:length(faxNoPadd),ixEpoch))));
plot(faxNoPadd,20*log10((XenopadR(1:length(faxNoPadd)))));
ylabel('Power [dB]')
xlabel('Frequency [Hz]')
linkaxes(ax5,'x')

%% Reconstruct the time waveform
vFixed=gs_istft(sdftrFixed,framesFixed);
vEpoch=gs_istft(sdftrEpoch,framesEpoch);
vr=gs_istft_nzp(stftr,framesEpoch);


figure(9); clf;
%stem(framesEpoch/fs,max(s)*ones(size(framesEpoch))); hold on
ax9(1)=subplot(311);
plot(tt,s,'k');hold on;
plot(tt(1:length(vFixed)),vFixed);
ax9(2)=subplot(312);
plot(tt,s,'k');hold on;
plot(tt(1:length(vEpoch)),vEpoch);
ax9(3)=subplot(313);
plot(tt,s,'k');hold on;
plot(tt(1:length(vr)),vr);
linkaxes(ax9,'x')
xlim([1.49 1.54]);
xlabel('Time [s]')
%ylabel('Normalised speech signal amplitude');

%% Calculate SNR
SNRfixed=calculateSNR(s,vFixed);
SNRepoch=calculateSNR(s,vEpoch);
SNRepochNoPad=calculateSNR(s,vr);

disp([SNRfixed, SNRepoch, SNRepochNoPad]);




