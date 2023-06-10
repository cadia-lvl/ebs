par=projParam();

%load RawRes.mat
%out1=out{1};
snrCube=zeros(length(par.db.numMels),length(par.db.frameMethod),length(out));

ii=1;
[nrw,ncl]=size(out);
for in=1:nrw
    for im=1:ncl
        outDr=out{in,im};
        for is=1:length(outDr)
           outUtt=outDr{is};
           snrCube(:,:,ii)=outUtt.snrMat;
           pesqCube(:,:,ii)=outUtt.pesqMat;
           outSeq{ii}=outUtt;
           ii=ii+1;
        end;
    end;
end
m_snrMat=mean(snrCube,3);
s_snrMat=std(snrCube,0,3);
m_pesqMat=mean(pesqCube,3);
s_pesqMat=std(pesqCube,0,3);

%% 
LabelFontSize=24;
TickFontSize=16;

figure(1); clf;
ax1(1)=subplot(2,1,1);
for ir=1:3
    errorbar(par.db.numMels,m_snrMat(:,ir),s_snrMat(:,ir),'LineWidth',1); hold on;
end;
grid on;
ylabel('SNR [dB]','FontSize',LabelFontSize);
set(gca,'FontSize',TickFontSize)
legend({'F','E','A'},'Location','NorthWest');
ax1(2)=subplot(2,1,2);
for ir=1:3
    errorbar(par.db.numMels,m_pesqMat(:,ir),s_pesqMat(:,ir),'LineWidth',1); hold on;
end;
grid on;
ylabel('PESQ score','FontSize',LabelFontSize);
xlabel('Number of Mel filters','FontSize',LabelFontSize);
set(gca,'FontSize',TickFontSize)

