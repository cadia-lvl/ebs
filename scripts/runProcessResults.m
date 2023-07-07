par=projParam();

%load RawRes.mat
%out1=out{1};
snrCube=zeros(length(par.db.numMels),length(par.db.config),length(out));

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

rgroup(1).config=[1 2 3 4];
rgroup(1).labels={'F','Fne','E','A'};
rgroup(1).title='No DC correction, with Zero Padding';

rgroup(2).config=[1 5 6 7] ;
rgroup(2).labels={'F','Fne','E','A'};
rgroup(2).title='With DC correction, with Zero Padding';

rgroup(3).config=[1 8 9 10];
rgroup(3).labels={'F','Fne','E','A'};
rgroup(3).title='With DC correction, no Zero Padding';

for irr=1:length(rgroup)

    figure(irr); clf;
    ax1(1)=subplot(2,1,1);
    for is=1:length(rgroup(irr).config)
        ir=rgroup(irr).config(is);
        errorbar(par.db.numMels,m_snrMat(:,ir),s_snrMat(:,ir),'LineWidth',1); hold on;
    end;
    grid on;
    title(rgroup(irr).title);
    ylabel('SNR [dB]','FontSize',LabelFontSize);
    set(gca,'FontSize',TickFontSize)
    legend(rgroup(irr).labels,'Location','NorthWest');
    ax1(2)=subplot(2,1,2);
    for is=1:length(rgroup(irr).config)
        ir=rgroup(irr).config(is);
        errorbar(par.db.numMels,m_pesqMat(:,ir),s_pesqMat(:,ir),'LineWidth',1); hold on;
    end;
    grid on;
    ylabel('PESQ score','FontSize',LabelFontSize);
    xlabel('Number of Mel filters','FontSize',LabelFontSize);
    set(gca,'FontSize',TickFontSize)
end

