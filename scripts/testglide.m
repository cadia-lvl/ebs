% testglide: glide example with constant frequency section
%%
addpath ../utils/
addpath ../tools/

% Parameters
%
fs=16000; % sample frequency
ts=0.5;                     % speech signal length in seconds
fxx=[100 190];               % range of fundamental frequencies
aww=[6 25]*1e-3;            % spectrogram analysis window lengths
fmax=3000; % frquency range of spectrogram
xlim=[0 0.5]; % time range to plot
rampstart=0.2; % time at which to start the frequency ramp
%
par=projParam;              % define computer-dependent parameters
pardb=par.db;
maxoff=pardb.maxoff;        % maximum allowed offset
%
% generate signal
%
load('lpcar');              % load the LPC parameters: ar=LPC filter, tf=timit file name, tlpc(1,2)=interval for LPC
ugt=(xlim(1)-0.1:1/fs:xlim(2)+0.1)';           % time vector extends 0.1 beyond xlim at each end
ns=length(ugt);             % number of speech samples
ugf=fxx(1)+(fxx(2)-fxx(1))*max(min((ugt-rampstart)/(xlim(2)-rampstart),1),0); % instantaneous pitch
ugph=cumsum(ugf/fs);        % cumulative phase
ug=v_glotlf(0,ugph);        % glottal flow waveform
s=filter([1 -0.95],ar,ug.*ugf.^(-1));  % apply ar filter plus lip radiation differentiator plus 1/fx gain correction
%
% now determine the frame limits
%
display(pardb.GCImethod);
periodlimk=round(fs./pardb.pitchlim);                     % {max targt mmin} pitch periods (samples)
framek=gs_frames(s,fs,pardb); % find frame boundaries
[framex,costgaindb,finaloff]=gs_frameadj(framek,s,fs,maxoff);
framekk=[[1 framex(1:end-1)+1];framex];             %  start and end samples for each frame
if any(framekk(:)<1) || any(framekk(:)>ns)
    error('frame boundaries do not lie within 1:ns');
end
framelen=framekk(2,:)-framekk(1,:)+1;
if (any(framelen>periodlimk(1)))
    error('Some frames are too long')
end
framemid=[0.5 0.5]*framekk/fs;                      % times of frame centres in seconds
%
% now perform transform
%
nfft=2*ceil(periodlimk(1)/2);                       % length of fft (force to be an even number)
[sdft,fax]=gs_stft(s,framex,nfft);
fax=fax*fs;                                         % convert frequency axis to Hz
nfftp=length(fax);                                  % number of positive frequencies
nframe=size(sdft,2);
%
% interpolate power spectrum onto regular grid (by replication)
%
sps=abs(sdft(1:nfftp,:)').^2; % single-sided power spectrum
staxs=(1:periodlimk(3):ns)'; % uniform time axis in samples
nstax=length(staxs);
[dum,ix]=sort([staxs; framex(:)]);  % compare with frame end samples
jx=1:length(ix);
jx(ix)=jx-jx(ix)+1;                 % sdft frame for each value in staxs
jx=min(jx(1:nstax),nframe);         % force <= nframe just in case
spsu=sps(jx,:);                     % power spectrum replicated onto uniform time axis
spsu(:,1)=0;                        % force DC component to zero

%%
% Now do plots
LabelFontSize=24;
TickFontSize=20;
PosVector=[ 36   425   985   522];
lwidth=1.2;

%
% Plot time waveform + frame boundaries
%
figure(2);
plot(ugt,s,'-b',ugt([1 end]),[0 0],':k');
hold on
v_axisenlarge([-1 -1.05]);
ylim=get(gca,'ylim');
plot(([framex;framex]-0.5)/fs+ugt(1),ylim,':k',(framek-1)/fs+ugt(1),0,'xr'); % plot frame boundaries (midway between samples)
hold off
xlabel('Time (s)');
ylabel('s(t), Unadj Frames (x), Frames');
%
% Plot instantaneous frequency
%
figure(4);
plot(ugt,ugf);
v_axisenlarge([-1 -1.05]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)');
% figure(3);
% fmaxi=floor(1+fmax/fax(2));           % number of frequency bins <= fmax
% spgrambw(spsu(:,1:fmaxi),[fs/periodlimk(3) ugt(1) fax(2) 0],'pJci');
% for i=1:length(aww)
% figure(10+i);
% spgrambw(s,fs,'pJcwi',1.81/aww(i),fmax);
% title(sprintf('Spectrogram: analysis window = %.1f ms',1000*aww(i)));
% end
naww=length(aww);
figure(1); clf; set(gcf,'Position',[90   241   656   706])
subplot(3+2*naww,1,1);
plot(ugt,s,'-k');
axis off
subplot(3+2*naww,1,3+2*naww+[-1 0]);
fmaxi=floor(1+fmax/fax(2));           % number of frequency bins <= fmax
v_spgrambw(spsu(:,1:fmaxi),[fs/periodlimk(3) ugt(1) fax(2) 0],'pJi');
set(gca,'xlim',xlim,'xtick',0:0.1:0.5,'xticklabel',{'0','0.1','0.2','0.3','0.4','0.5'});
v_texthvc(0.03,0.95,'Epoch-based frames','LTk');
ylabel('kHz');
chca=get(gca,'Children');
set(chca(1),'FontSize',TickFontSize);
xlabel('Time [s]');
yticks(1000*[0, 1, 2, 3]); yticklabels([0, 1, 2, 3]);
set(gca,'FontSize',TickFontSize)
for i=1:length(aww)
subplot(3+2*naww,1,2*i:2*i+1);
v_spgrambw(s,[fs ugt(1)],'pJi',1.81/aww(i),fmax);
yticks(1000*[0, 1, 2, 3]); yticklabels([0, 1, 2, 3]);
set(gca,'xlim',xlim,'xtick',0:0.1:0.5,'xticklabel',[],'xlabel',[]);
ylabel('kHz');
v_texthvc(0.03,0.95,sprintf('Fixed %.0f ms frames',1000*aww(i)),'LTk');
chca=get(gca,'Children');
set(chca(1),'FontSize',TickFontSize);
set(gca,'FontSize',TickFontSize)
end
%print('./Figures/spectroGrams.pdf','-dpdf');