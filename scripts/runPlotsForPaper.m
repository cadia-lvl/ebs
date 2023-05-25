%Plot parameters
LabelFontSize=24;
TickFontSize=20;
PosVector=[ 36   425   985   522];
lwidth=1.2;


par=projParam();
spfile=[par.pth.speechpth '/timit/test/dr3/mjvw0/sx113.wav'];

[sp,fs]=v_readsph(spfile,'wt');
tt=(0:length(sp)-1)/fs;

phn=extractPhn([spfile(1:end-3) 'PHN']);
[framek,gci]=gs_frames(sp,fs);

figure(1); clf; set(gcf,'Position', PosVector);
subplot(3,1,1);
plot(1000*tt,sp,'k','LineWidth',lwidth); hold on;
stem(1000*gci,0.04*ones(size(gci)),'b','LineWidth',lwidth)
xlim([550  730])
ylabel('Voiced','FontSize',LabelFontSize)
set(gca,'YTick',[])
set(gca,'FontSize',TickFontSize)


subplot(3,1,2);
plot(1000*tt,sp,'k','LineWidth',lwidth); hold on;
stem(1000*gci,0.02*ones(size(gci)),'b','LineWidth',lwidth)
xlim([1736    2008])
ylabel('Unvoiced','FontSize',LabelFontSize)
set(gca,'YTick',[])
set(gca,'FontSize',TickFontSize)


subplot(3,1,3);
plot(1000*tt,sp,'k','LineWidth',lwidth); hold on;
stem(1000*gci,4e-3*ones(size(gci)),'b','LineWidth',lwidth)
xlim([1158   1360])
ylabel('Short pause','FontSize',LabelFontSize)
xlabel('Time [ms]','FontSize',LabelFontSize)
set(gca,'YTick',[])
set(gca,'FontSize',TickFontSize)
set(gcf,'PaperOrientation','landscape');


%%

framesF=getFrames(sp,fs,par.db,1);  % Fixed frames
framesE=getFrames(sp,fs,par.db,2);  % Epoch frames
framesM=getFrames(sp,fs,par.db,3);  % Epoch Modified frames

% Picking out
ixF=find(698<1000*framesF/fs & 1000*framesF/fs <705);
ixM=find(704<1000*framesM/fs & 1000*framesM/fs <705);
ixE=find(704<1000*framesE/fs & 1000*framesE/fs <705);


[sdftF,faxF]=gs_stft(sp,framesF,par.db.nfft);
[sdftE,faxE]=gs_stft(sp,framesE,par.db.nfft);
[sdftM,faxM]=gs_stft(sp,framesM,par.db.nfft);

frame=round((699.188/1000)*fs):(round((721.938/1000)*fs)-1);
sframe=sp(frame);
Sframe=fft(sframe);
ff=fs*(0:length(Sframe)-1)/length(Sframe);
ff=ff(1:floor(length(Sframe)/2));
Sframe=Sframe(1:floor(length(Sframe)/2));



%%

figure(2); clf; set(gcf,'Position', [36   299   733   648]);
subplot(3,1,1);
stem(1000*framesF/fs,1.2*ones(size(framesF))*max(sp),'sq--','LineWidth',lwidth); hold on;
stem(1000*framesE/fs,0.7*ones(size(framesE))*max(sp),'LineWidth',lwidth)
stem(1000*framesM/fs,-0.7*ones(size(framesM))*max(sp),'LineWidth',lwidth)
stem(1000*framesF/fs,-ones(size(framesF))*max(sp),'bsq--','LineWidth',lwidth)

plot(1000*tt,sp,'k','LineWidth',lwidth);

ylabel('Speech','FontSize',LabelFontSize)
xlabel('Time [ms]','FontSize',LabelFontSize)
set(gca,'YTick',[])
set(gca,'FontSize',TickFontSize)
axis([695  730  -0.06    0.08])

subplot(3,1,2);
plot(fs*faxF/1000,20*log10(abs(sdftF(1:length(faxF),ixF))),'LineWidth',lwidth); hold on;
plot(fs*faxE/1000,20*log10(abs(sdftE(1:length(faxE),ixE))),'LineWidth',lwidth) 
plot(fs*faxM/1000,20*log10(abs(sdftM(1:length(faxM),ixM))),'LineWidth',lwidth)

ylabel('Power [dB]','FontSize',LabelFontSize)
grid on;
set(gca,'FontSize',TickFontSize)
axis([0  5  -50    20])

subplot(3,1,3);
plot(fs*faxF/1000,unwrap(angle(sdftF(1:length(faxF),ixF))),'LineWidth',lwidth); hold on;
plot(fs*faxE/1000,unwrap(angle(sdftE(1:length(faxE),ixE))),'LineWidth',lwidth) 
plot(fs*faxM/1000,unwrap(angle(sdftM(1:length(faxM),ixM))),'LineWidth',lwidth)

ylabel('Phase [rad]','FontSize',LabelFontSize)
xlabel('Frequency [kHz]','FontSize',LabelFontSize)
grid on;
set(gca,'FontSize',TickFontSize)
set(gcf,'PaperOrientation','landscape');
axis([0  5  -100    15])

