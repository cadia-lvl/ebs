%Plot parameters
LabelFontSize=24;
TickFontSize=20;
PosVector=[56          71        1107         884];

%% Importing file and extracting GCIs
%
addpath VoiceSourceTools
%spfile='jg_ae_fm_vm_sp';  

% From ColdFish - select between 3.5s and 4.1 male voice saying "Eg man"
spfile='fromColdFish';
%spfile='fromColdFish';
[sp,fs] = readwav(spfile); tt=(0:length(sp)-1)/fs; nn=3.5<tt&tt<4.1; sp=sp(nn);

tt=(0:length(sp)-1)/fs;
[gci, goi] = dypsagoi(sp,fs);

% choose smaller segment
tb=0.36;
te=tb+0.08;
tc=tb+0.04;

%% Copy-paste a cycle for analysis

gFrame=1000*tt(gci(42:43));  %Example glottal cycle
sp1=sp(gci(42):(gci(43)-1));

spSyn4=repmat(sp1,4,1); ttSyn4=(0:length(spSyn4)-1)/fs;
spSyn8=repmat(sp1,8,1); ttSyn8=(0:length(spSyn8)-1)/fs;
spSyn11=repmat(sp1,10,1); ttSyn11=(0:length(spSyn11)-1)/fs;


N1=length(sp1);
a1=fft(sp1)/N1;
ffgseg=(fs)*(0:(N1-1))/N1;

N4=length(spSyn4);
a4=fft(spSyn4)/N4;
ffsyn4=(fs)*(0:(N4-1))/N4;

N8=length(spSyn8);
a8=fft(spSyn8)/N8;
ffsyn8=(fs)*(0:(N8-1))/N8;

aa1=abs(a1);
aa4=abs(a4);
aa8=abs(a8);

laa1=10*log10(abs(a1));
laa4=10*log10(abs(a4));
laa8=10*log10(abs(a8));


%% Exact Gcycle analysis
sp4=sp(gci(42):(gci(46)-1));
sp8=sp(gci(42):(gci(50)-1));

N4r=length(sp4);
a4r=fft(sp4)/N4r;
ffr4=(fs)*(0:(N4r-1))/N4r;

N8r=length(sp8);
a8r=fft(sp8)/N8r;
ffr8=(fs)*(0:(N8r-1))/N8r;

aa4r=abs(a4r);
aa8r=abs(a8r);

laa4r=10*log10(abs(a4r));
laa8r=10*log10(abs(a8r));

% Get average (and std) glottal cycle durations in plot:
dg=diff(gci(41:51));
mdg=1000*mean(dg)/fs;
sdg=1000*std(dg)/fs;

%% Resample Gcycle analysis

Nseg=N1;
sp4s=zeros(4*Nseg,1);
sp4s(1:Nseg)=sp1;
for ik=1:4
    spseg=sp(gci(42+ik):(gci(42+ik+1)-1));
    sp4s((1:Nseg)+ik*Nseg)=resample(spseg,Nseg,length(spseg));  
end
sp8s=zeros(4*Nseg,1);
sp8s(1:Nseg)=sp1;
for ik=1:4
    spseg=sp(gci(42+ik):(gci(42+ik+1)-1));
    sp8s((1:Nseg)+ik*Nseg)=resample(spseg,Nseg,length(spseg));  
end

N4s=length(sp4s);
a4s=fft(sp4s)/N4s;
ffs4=(fs)*(0:(N4s-1))/N4s;

N8s=length(sp8s);
a8s=fft(sp8s)/N8s;
ffs8=(fs)*(0:(N8s-1))/N8s;

aa4s=abs(a4s);
aa8s=abs(a8s);

laa4s=10*log10(abs(a4s));
laa8s=10*log10(abs(a8s));

%% "Normal" frame - spectrum of
nfft=1024;
wsiz=round(25e-3*fs);
hsiz=round(10e-3*fs);

[Sk,ff,tT] = spectrogram(sp,hamming(wsiz),wsiz-hsiz,nfft,fs,'yaxis');
tIx=39;
tFrame=1000*[tT(tIx)-wsiz/(2*fs) tT(tIx)+wsiz/(2*fs)];
ntFrame=tFrame*fs/1000;



%% Two panels of time traces

f1=figure(1); clf; set(gcf,'Position',PosVector);
subplot(311);
plot(1000*tt,sp,'k'); hold on;
%stem(1000*tt(gci),max(sp)*ones(size(gci))); yl1=ylim;
yl1=[-1 1];
patch(1000*[tb te te tb],[yl1(1) yl1(1) yl1(2) yl1(2)],'k','FaceAlpha',0.1,'Edgecolor','none')
yticks([])
ylabel('Speech signal','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)

nn=find(tb<tt&tt<te);
subplot(312);
%plot(1000*tt(nn),sp(nn),'k'); hold on;
plot(1000*tt,sp,'k'); hold on;
stem(1000*tt(gci),ones(size(gci)),'Color',0.7*[1 1 1],'LineWidth',1.5,'Marker','none','BaseValue',-1); hold on;
yl1=[-1 1];
patch([tFrame tFrame(end:-1:1)],[yl1(1) yl1(1) yl1(2) yl1(2)],'k','FaceAlpha',0.2,'Edgecolor','none'); hold on;
patch([gFrame gFrame(end:-1:1)],[yl1(1) yl1(1) yl1(2) yl1(2)],'k','FaceAlpha',0.05,'Edgecolor','none')
yticks([])
xlim(1000*[tb te])
ylabel('Speech signal','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)

subplot(313);
ttc=1000*(0:1023)/fs;
plot(ttc(1:wsiz),hamming(wsiz).*sp(ntFrame(1):(ntFrame(2)-1)),'k'); hold on;
plot(ttc(1:wsiz),hamming(wsiz),'k--')
plot(ttc(wsiz+1:1024), zeros(1,524),'k')
yticks([])
xlim(1000*[0 1023]/fs)
ylabel('Analysis signal','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)
xlabel('Time [ms]','FontSize',LabelFontSize);

f11=figure(11); clf; set(gcf,'Position',PosVector);
subplot(313);
%plot(1000*tt(nn),sp(nn),'k'); hold on;
plot(1000*ttSyn11,spSyn11,'k'); hold on;
%stem(1000*tt(gci),ones(size(gci)),'Color',0.7*[1 1 1],'LineWidth',1.5,'Marker','none','BaseValue',-1); hold on;
yl1=[-1 1]; ylim(yl1);
xlim(1000*[0 ttSyn11(end)])
yticks([])
ylabel('Replicated cycle','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)

xlabel('Time [ms]','FontSize',LabelFontSize);


plotName='Figures/signal1';
hgsave(f1, [plotName '.fig']);
print(f1, [plotName '.eps'],'-deps');


%% Spectrum of

f2=figure(2); clf; set(gcf,'Position',PosVector);
subplot(311);
plot(ffsyn8/1000,aa8,'-','Color',0.4*[1 1 1],'LineWidth',1.5); hold on;
plot(ffsyn4/1000,aa4,'sq--','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
plot(ffgseg/1000,aa1,'o:','Color',0*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
ylabel('Repeated cycle','FontSize',LabelFontSize)
%ylabel('Power/freq [dB/Hz]','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)
xlim([0 1]);
subplot(312);
plot(ffs8/1000,aa8s,'-','Color',0.4*[1 1 1],'LineWidth',1.5); hold on;
plot(ffs4/1000,aa4s,'sq--','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
plot(ffgseg/1000,aa1,'o:','Color',0*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
ylabel('Cycles resampled','FontSize',LabelFontSize)
%ylabel('Power/freq [dB/Hz]','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)
xlim([0 1])
subplot(313);
plot(ffr8/1000,aa8r,'-','Color',0.4*[1 1 1],'LineWidth',1.5); hold on;
plot(ffr4/1000,aa4r,'sq--','Color',0.3*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
plot(ffgseg/1000,aa1,'o:','Color',0*[1 1 1],'LineWidth',1.5,'MarkerSize',10);
ylabel('Different cycles','FontSize',LabelFontSize)
%ylabel('Power/freq [dB/Hz]','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)
xlim([0 1])
xlabel('Frequency [kHz]','FontSize',LabelFontSize)

plotName='Figures/aks';
hgsave(f2, [plotName '.fig']);
print(f2, [plotName '.eps'],'-deps');




%% Spectrogram using a "normal" frame

fpor=4;
f3=figure(3); clf; set(gcf,'Position',PosVector);
subplot(212);
surf(1000*tT,ff(1:(nfft/fpor))/1000,10*log10(abs(Sk(1:(nfft/fpor),:))),'EdgeColor','none'); view(0,90), axis xy; axis tight; colormap((bone)); 
hold on;
plot(1000*tT(tIx)*[1 1],[0 fs/fpor]/1000,'w');
plot(1000*tT(tIx+1)*[1 1],[0 fs/fpor]/1000,'w');
ylabel('Frequency [kHz]','FontSize',LabelFontSize)
xlabel('Time [ms]','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)
subplot(211);
plot(ff(1:(nfft/fpor))/1000,10*log10(abs(Sk(1:(nfft/fpor),tIx))),'k');
grid on
%yticks([])
ylabel('Power/freq [dB/Hz]','FontSize',LabelFontSize)
xlabel('Frequency [kHz]','FontSize',LabelFontSize)
set(gca,'FontSize',TickFontSize)

plotName='Figures/spectrogram1';
hgsave(f3, [plotName '.fig']);
print(f3, [plotName '.eps'],'-deps');

 



