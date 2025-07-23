% clear;
close all;
timit=gettimitpath;                         % get path to timit subfolder of timit database
figsuf=char(string(datetime("now","Format","yyy-MM-dd_hh-mm")));                % current date and time as suffix for output filenames
figpdf=['figures_<m>/<m>_<n>_' figsuf];
% newdata='TEST/DR3/MMAB0/SX282.WAV';       % 'ep pwlin-ph' gives PESQ<4.2 @ 200 mel bins
% newdata='TEST/DR4/FCRH0/SX278.WAV';       % /s/ @ t=3.1s gives strange framelengths with par.smoothalg='quadlog' because of a spurious frame @ t=3.04
newdata='TRAIN/DR6/MSJK0/SX246.WAV';        % demo file for documentation
% newdata=0;                                % create new data for each run (1 or 0 or tf output from previous run)
% algorithm parameters
parvar={'cfname'};                          % cell array row listing the position-dependent parameters within configpars structure
configpars={                                % one row per trial giving the parameter values
    {'30ms'     'epoch' 0 'fixfl'  0.03};
    {'6ms'     'epoch' 0 'fixfl'  0.006};
    {'ep+q' 'smoothalg' 'quadlog'};
    {'ep+q-natc' 'smoothalg' 'quadlog' 'interpstft' 'natural'};
    {'ep+q-natm' 'smoothalg' 'quadlog' 'interpstft' 'natural' 'interpdom' 'magcph'};
    };
nparvar=length(parvar); % number of position-dependent parameters
nconfig=length(configpars);                   % number of configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plots (include in plotlist to plot):
%           1       Spectrogram
%           2       Configurations
%           3       Frame lengths
%           4       Phone labels, durations and phone-type classification
%           5       Plot of magnitude change between consecutive frames
%           6       Plot of complex-value change between consecutive frames
%           7       Plot of complex-value change between consecutive frames interpolating in magnitude and unwrapped-phase
%           8       Plot of benefit of interpolating in magnitude and unwrapped-phase; difference between plots 6 and 5
%           9       Horizontal slice through "spectrogram" (fig 100) at frequency fplot
%          10       Vertical slice through "spectrogram" (fig 100) at time tplot
%          11       Bar-graph of complex-value change between consecutive frames interpolating in magnitude and unwrapped-phase for each phone-type
%          12       Bar-graph of magnitude change between consecutive frames for each phone-type
%          13       Bar-graph of complex-value change between consecutive frames for each phone-type
%          14       Bar-graph of within-phone magnitude variance
%          15       Bar-graph of within-phone complex variance
%    ... plots for each configuration (added to configuration number).
%         100+      "spectrogram" of STFT
%         200+      "spectrogram" of interpolated STFT (only if par.interpstft~='none')
%         300+      Bhattacharyya divergence within each phone
%
plotlist=[1:4 9:10 12:15 100];       % list of graphs to plot
dopdf=0;                            % set to 1 to export pdfs of figures
figbold=0;                          % set to 1 to embolden selected figures
playconfig=[];                      % list of configurations to play; 0=original
tplot=0.28;                         % interesting time to plot [seconds]
fplot=1000;                         % interesting frequency to plot [Hz]
%%% mapping of phone classes (from w_phoncode.m) to a small number of phone types
phncls2typ=[6 6 4 4 4 4 4 4 4 4 4 4 4 2 3 4 4 4 6 6 6 6 6 6 6 6]; % [1=changed] 2=vowel, 3=dipthong, 4=consonant, [5=unvoiced], 6=other
phntypename={'change','vowel', 'dipth', 'v-cons', 'u-cons', 'other'};
nphntype=length(phntypename);       % number of phone types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Default parameter values     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% = global =
par0.epoch=       1;                 % 0=fixed frames, 1=epoch-based frames, 2=epoch-based interpolated in log-domain to fixed
par0.fixfl=       0.03;              % target frame length in seconds [0.03]
par0.autoscale=   'energy';          % scale reconstructed waveform to match: {'none','energy'}
par0.cfname=      '';                % configuration name
% = gs_frames =
par0.pitchlim=    [40 50 400];       % [min targt max] pitch (Hz)
par0.gcifrac=     0.3;               % position of GCI in analysis frame
par0.GCImethod=   'YAGA';            % GCI estimation method
par0.ovfact=      1;                 % integer overlap factor; 1 for no overlap, 2 for 50%
% = smoothframes
par0.smoothalg=   'none';            % smoothing algorithm {'none','lin','quadlin','quadlog'}
par0.voithresh=   0.7;               % cross-correlation threshold for voiced frames
par0.voilenrat=   1.1;               % length ratio threshold for voiced frames
% = stfte =
par0.window=      'r';               % Window: {'r','n','m'} = {Rectangular, Hanning, Hamming}
par0.offset=      'none';            % offset removal: {'none','mean'}
par0.scale=       'none';            % scaling method: {'none','peakabs','rms'}
par0.pad=         'none';            % zero-padding method: {'none','zero','ends'}
par0.groupdelay=  'none';            % linear phase component: {'none','dct','ewgd','ewgdint','cplx','cplxint','xcor'}
% = stftegrid =
par0.interpstft=  'none';            % interpolation method for call to griddata: {'none','indep','nearest','linear','natural','cubic','v4'}
par0.interpflen=  0.007;             % target frame length in seconds (reciprocal of frequency-grid spacing) [0.01]
par0.interpov=    1;                 % overlap factor; 1 for no overlap, 2 for 50%; time grid-spacing is par.interpflen/par.interpov
par0.interpdom=   'cplx';            % interpolation domain: {'cplx','magcph','crmcph'}
par0.interpfsps=   1e-5;             % interpolation scale factor: fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second (fs per sample)
par0.interpext=    'omit';           % Handling of extrapolated frames: {'omit','zero','rep','refl'}
% = stfte2melc
par0.fbank=       'm';               % filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel}
par0.keepDC=      1;                 % preserve DC as the lowest MEL bin {0, 1}
par0.MELphase=    'true';            % MEL STFT phase calculation: {'true','zero','linear','piecewiselin','grpdel'}
par0.MELdom=      'mag';             % MEL filterbank domain: {'mag', 'pow'}
par0.regwt=       0.01;              % regularizing factor for phase weights when par.MELphase='piecewiselin'
par0.loops=       5;                 % maximum number of iterations when par.MELphase='piecewiselin'
par0.ssqthr=      1.0;               % stopping threshold when par.MELphase='piecewiselin'
% = melc2stfte
par0.invmethod=  'pinv';             %  MEL inversion method: {'pinv','transp','linterp'} ['linterp']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read speech data             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(newdata)                                              % if newdata is a character string ...
    tf=newdata;                                                 % it specifies a file path
    [s,fs,wrd,phn]=gettimit(tf,'n',timit);                      % read the file
elseif newdata==1 || ~exist('s','var')                          % check if new data needed
    [tf,ty,tk,s,fs,wrd,phn]=timitfiles('pz',1,timit);           % read a level-normalized random TIMIT file
end
tplot=min(tplot,length(s)/fs);                                  % ensure that plot time lies within the signal
nphn=size(phn,1);                                               % number of phones in sentence
for i=1:nphn
    [phn{i,3},phn{i,4}]=w_phoncode('tU',phn{i,2});              % append unicode and distinctive feature information
end
phnlim=cell2mat(phn(:,1));                                      % extract phone limits (in seconds)
phnfeat=cell2mat(phn(:,4));                                     % extract phone distinctive features
phntyp=phncls2typ(phnfeat(:,3)); % map phone class to one of 5 types
phntyp(bitand(phnfeat(:,5),2^11)==0 & (phntyp(:)==4))=5; % set unvoiced consonants to 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Initialize results array     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfnames=cell(nconfig,1);                                    % list of configuration names
cfresults=cell(nconfig,8);                                  % space for meta information
framephntcnts=zeros(nphntype,3,nconfig);                    % space for phone counts
axlinkt=[];                                                 % time axis linking list
axlinkf=[];                                                 % frequency axis linking list
axlinknm=[];                                                % number of mel bins linking list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Loop through each configuration   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icfg=1:nconfig                                          % loop through parameter configurations
    [par,configtxt]=parupdate(par0,configpars{icfg},parvar,'%cfname%: $c&=%&%$, $');
    cfnames{icfg}=par.cfname;                                   % save configuration name
    nhopf=round(par.interpflen*fs/par.interpov);                % frame hop for fixed frames in samples
    nbinf=2*round(par.interpflen*fs/2);                         % effective dft length (always even) for fixed frames
    interphzps=par.interpfsps*fs^2;                             % convert interpolation scale factor from fs/sample to Hz/s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  process the current configuration  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_finishat([icfg 1 nconfig]); % one row per nested loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        perform framing              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [frendraw,gci]=gs_frames(s,fs,par);                         % index of last sample in each frame
    gcik=round(gci*fs+1);                                       % gci positions (sample numbers)
    if par.epoch                                                % if epoch-based frames
        frameend=smoothframes(s,frendraw,par);                  % smooth the frame pitch contour
    else                                                        % if fixed length frames
        hop=round(par.fixfl*fs/par.ovfact);                     % force frame hop to be an integer number of samples
        frameend=hop:hop:length(s);                             % ends of fixed frames of length hop
    end
    framelim=[[1 frameend(1:end-par.ovfact)+1];frameend(par.ovfact:end)];   % start and end of each frame (sample indices)
    metain=[framelim(1,:);1+[-1 1]*framelim]';                  % one row per frame: [start-sample frame-length]
    [stft,meta,grpd]=stfte(s,metain,[],par);                    % epoch-based STFT
    % [stft,meta]=stftgrid(stft,meta,par);                      % OLD: optionally map onto a fixed grid and update stft and meta
    grid=[meta(1,1) round(par.interpflen*fs) round(par.interpflen*fs) round(par.interpflen*fs/par.interpov)];      % [firstsamp framelen ndft nhop] for uniform grid
    [stft,meta]=stftegrid(stft,meta,grid,par);                  % optionally map onto a fixed grid
    tax=(meta(:,1:2)*[1;0.5]-0.5)/fs;                           % time axis (frame centres in seconds, 1st sample @ 1/fs)
    nfftmax=size(stft,2);
    nfftp=1+floor(nfftmax/2); % number of positive frequencies
    fax=(0:nfftp-1)*fs/nfftmax; % frequency bins in Hz
    framephn=v_interval(tax(:),[phnlim(:,1);phnlim(end,2)],'cC');  % determine phone index for each frame
    framephnlim=cumsum(full(sparse(framephn,1,1)));
    framephnlim=[[1;framephnlim(1:end-1)+1] framephnlim]; % first and last frames corresponding to each phone (empty if last<first)
    framephntype=phntyp(framephn);   % map to a reduced set of phone types (as listed in phncls2typ)
    framephnspurt=[true framephntype(2:end)~=framephntype(1:end-1)]; % flag the start of each phone-type spurt
    framephntcnts(:,1:2,icfg)=[full(sparse(framephntype(framephnspurt),1,1,nphntype,1)) full(sparse(framephntype,1,1,nphntype,1))]; % count of spurts and frames by type
    framephntcnts(:,3,icfg)=framephntcnts(:,2,icfg)-framephntcnts(:,1,icfg);
    framephntcnts(1,3,icfg)=sum(framephntcnts(:,1,icfg))-1;               % framephntcnts cols: 1=# spurts, 2=# frames, 3=# frame-pairs
    nframe=size(meta,1);                                    % number of frames
    [dum,tploti]=min(abs(tax-tplot)); % index of frame to plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   calculate inter-frame increments  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delstft=zeros(nframe-1,3); % inter-frame increments in dB
    dbclip=-20; % clip errors to >=-20 dB
    for jframe=2:nframe
        iframe=jframe-1; % base frame
        stfti=stft(iframe,1:meta(iframe,2)); % extract frame i from stft
        stftia=abs(stfti);
        stftiu=unwrap(angle(stfti)); % unwrapped phase
        eni=sum(stftia.^2); % energy of frame i
        stftj=stft(jframe,1:meta(jframe,2)); % extract frame j from stft
        stftja=abs(stftj);
        stftju=unwrap(angle(stftj));
        if meta(jframe,2)==meta(iframe,2) % consecutive frames are the same length
            delstft(iframe,1)=max(db(sum((stftja-stftia).^2)/eni)/2,dbclip);
            delstft(iframe,2)=max(db(sum(abs(stftj-stfti).^2)/eni)/2,dbclip);
            delstft(iframe,3)=delstft(iframe,2); % same as complex since no interpolation
        else
            t=fftinterp(meta(jframe,2),meta(iframe,2))'; % sparse interpolation matrix
            delstft(iframe,1)=max(db(sum((stftja*t-stftia).^2)/eni)/2,dbclip); % interpolate magnitude only
            delstft(iframe,2)=max(db(sum(abs(stftj*t-stfti).^2)/eni)/2,dbclip);  % interpolate real-imag
            delstft(iframe,3)=max(db(sum(abs((stftja*t).*exp(1i*stftju*t)-stfti).^2)/eni)/2,dbclip);  % interpolate magnitude-phase
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate mean dB changes by phone type  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delstftphn=zeros(nphntype,size(delstft,2));
    for i=1:length(phntypename)
        if i==1
            msk=framephntype(1:end-1)~=framephntype(2:end);
        else
            msk=(framephntype(1:end-1)==i) & (framephntype(2:end)==i);
        end
        if any(msk)
            delstftphn(i,:)=mean(delstft(msk,:),1);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   calculate within-phone variances  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stftmpdb=zeros(nphn,2); % (:,1)=magnitude variances, (:,2)=complex variances
    if all(meta(:,3)==meta(1,3)) % if all ffts have the same length, so we can do averages
        for i=1:nphn
            framelist=framephnlim(i,1):framephnlim(i,2); % list of frames in this phone
            nframelist=length(framelist); % number of frames in this phone
            stftmav=mean(abs(stft(framelist,:)),1); % mean of magnitude spectrum within this phone
            stftmpdb(i,1)=mean(max(db(sum((abs(stft(framelist,:))-repmat(stftmav,nframelist,1)).^2,2)/sum(stftmav.^2))/2,dbclip)); % average dB error over phone
            stftcav=mean(stft(framelist,:),1); % mean of compex spectrum within this phone
            stftmpdb(i,2)=mean(max(db(sum(abs(stft(framelist,:)-repmat(stftcav,nframelist,1)).^2,2)/sum(abs(stftcav.^2)))/2,dbclip)); % average dB error over phone
        end
    end
    stftmpdbphn=zeros(nphntype,size(stftmpdb,2));
    for i=1:length(phntypename)
        msk=(phntyp'==i) & (framephnlim*[-1;1]>=0); % runs of this phone type excluding unsampled phones
        if any(msk)
            stftmpdbphn(i,:)=mean(stftmpdb(msk,:),1);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  save per-configuration information %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfresults{icfg,1}=meta;                                 % 1: metadata
    cfresults{icfg,2}=delstft;                              % 2: consecutive-frame differences in dB
    cfresults{icfg,3}=delstftphn;                           % 3: consecutive-frame differences in dB averaged over phone type
    tplotnf=meta(tploti,3); % DFT length for frame tploti
    tplotnfp=1+floor(tplotnf/2); % number of positive frequencies
    cfresults{icfg,5}=[(0:tplotnfp-1)'*fs/tplotnf stft(tploti,1:tplotnfp).'];      % 5: spectrum at single time: tplot
    cfresults{icfg,7}=stftmpdbphn;                          % 7: within-phone differences in dB averaged over phone type
    if all(meta(:,3)==meta(1,3))
        fploti=round(1+(fplot*nfftmax)/fs); % index of frequency to plot (fixed frames only)
        cfresults{icfg,4}=[tax(:) stft(:,fploti)];              % 4: trace of single frequency vs time: fplot
        cfresults{icfg,6}=[tax(tploti) fax(fploti) fax(fploti)];                      % 6: miscellaneous values
    else
        fploti=1+round(fplot/fs*meta(:,3)');
        fplotx=(fploti-1)*fs./meta(:,3)'; % actual frequency
        cfresults{icfg,4}=[tax(:) stft((1:nframe)+(fploti-1)*nframe).'];              % 4: trace of single frequency vs time: fplot
        cfresults{icfg,6}=[tax(tploti) min(fplotx) max(fplotx)];                      % 6: miscellaneous values
    end
    cfresults{icfg,8}=configtxt;                            % text description of configuration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Per-Configuration plots       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(plotlist==100)                                   % plot "spectrogram" of STFT
        figure(100+icfg);
        stfteplot(stft,meta,fs,40);
        if strcmp(par.interpstft,'none')
            title(cfnames{icfg});
        else
            title(sprintf('%s: %s@%.0fHz/hop',cfnames{icfg},par.interpstft,interphzps*nhopf/fs));
        end
        axlinkt=[axlinkt gca]; % link time axis
    end
    if any(plotlist==200) && ~strcmp(par.interpstft,'none') % plot spectrogram of interpolated STFT
        figure(200+icfg);
        imagesc(tax,fax,db(abs(stft(:,1:nfftp)))');
        axis 'xy';
        set(gca,'clim',get(gca,'clim')*[0;1]+[-40 0]);      % limit range to 40 dB
        hold on
        msk=faxv<=fs/2;
        plot(taxv,faxv,'xk');
        hold off
        colorbar;
        xlabel('Time (s)');
        ylabel('Freq (Hz)');
        if strcmp(par.interpstft,'none')
            title(cfnames{icfg});
        else
            title(sprintf('%s: %s@%.0fHz/hop',cfnames{icfg},par.interpstft,interphzps*nhopf/fs));
        end
        axlinkt=[axlinkt gca]; % link time axis
    end
    if any(plotlist==300) && ~isempty(dbfgm)
        figure(300+icfg);
        imagesc(dbfgm);
        colorbar;
        cblabel('Bhattacharyya div')
        set(gca,'clim',[0 30]);
        if strcmp(par.interpstft,'none')
            title(cfnames{icfg});
        else
            title(sprintf('%s: %s@%.0fHz/hop',cfnames{icfg},par.interpstft,interphzps*nhopf/fs));
        end
    end
end
delstftmu=reshape(cell2mat(cfresults(:,3)),nphntype,nconfig,3);
stftmpdbmu=reshape(cell2mat(cfresults(:,7)),nphntype,nconfig,2);
% play sound if required
if any(playconfig==0)
    soundsc(s,fs);
end
% global plots
for iplot=1:20
    cfnamesmu=cfnames;                          % initialize space for configuration names  extra plot-dependent information
    if any(plotlist==iplot)
        switch iplot
            case 1                              %%% plot 1: spectrogram %%%
                figure(iplot);
                v_spgrambw(s,fs,'pJcwaAtT',[],[],[],[],{wrd phn}); % plot spectrogram with word and phone transcriptions
                % hold on
                % plot(tplot*[1 1],get(gca,'ylim'),'-w');
                % hold off
                axlinkt=[axlinkt gca]; % link time axis
                title(['Original: ' tf]);
                if dopdf, v_fig2pdf(figpdf), end;

            case 2
                figure(iplot);
                cfdesc='';
                linelen=50;
                for i=1:nconfig
                    cfdesci=sprintf('%d) %s',i,cfresults{i,8});
                    while length(cfdesci)>linelen
                        j=find(cfdesci(5:linelen)==',',1,'last')+4;
                        if isempty(j)
                            j=find(cfdesci(5:linelen)=='=',1,'last')+4;
                            if isempty(j)
                                j=linelen;
                            end
                        end
                        cfdesc=[cfdesc cfdesci(1:j) char(10)]; % append remainder of description
                        cfdesci=['    ' cfdesci(j+1:end)];
                    end
                    if length(cfdesci)>0
                        cfdesc=[cfdesc cfdesci char(10)]; % append remainder of description
                    end
                end
                v_texthvc(0,0.95,cfdesc,'LTk');
                axis off
                title('Configurations');

            case 3   %%% plot 2: frame lengths %%%
                figure(iplot);
                subplot(4,1,2:4);
                for icfg=1:nconfig
                    meta=cfresults{icfg,1};
                    nframe=size(meta,1); % number of frames
                    plot(reshape([meta(:,1)';meta(:,1)'+meta(:,2)'-1],2*nframe,1)/fs,reshape(meta(:,[2 2])',2*nframe,1)*1000/fs);
                    hold on
                end
                hold off
                v_axisenlarge([-1 -1.05]);
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('Frame length (ms)');
                legend(cfnames,'location','best');
                subplot(4,1,1);
                plot((1:length(s))/fs,s);
                v_axisenlarge([-1 -1.05]);
                axlinkt=[axlinkt gca]; % link time axis
                title('Frame lengths');
                if dopdf, v_fig2pdf(figpdf), end;

            case 4 %%% 3: phone types
                figure(iplot);
                plot(reshape(phnlim',[],1),reshape(repmat(phntyp,2,1),[],1));
                hold on
                for i=1:size(phn,1)
                    v_texthvc(phnlim(i,:)*[0.5;0.5],phntyp(i),phn{i,3},'mbk');
                end
                for i=2:nphntype
                    v_texthvc(0.01,i,sprintf('%d',max(framephntcnts(i,1,:),[],3)),'Lmk');
                end
                hold off
                v_axisenlarge([-1 -1 -1.05 -1.1]);
                set(gca,'ytick',2:nphntype,'yticklabel',phntypename(2:nphntype));
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                title('Phone Types');
                if dopdf, v_fig2pdf(figpdf), end;

            case 5
                figure(iplot);
                for icfg=1:nconfig
                    meta=cfresults{icfg,1};
                    delstft=cfresults{icfg,2};
                    nframe=size(meta,1); % number of frames
                    plot((meta(1:nframe-1,1)'+meta(1:nframe-1,2)'-0.5)/fs,delstft(:,1));
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(delstft(:,1)))];
                    hold on
                end
                v_axisenlarge([-1 -1.05]);
                plot(get(gca,'xlim'),[0 0],':k');
                hold off
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('Inter-frame \Delta (dB)');
                legend(cfnamesmu,'location','best');
                title('Interframe Magnitude Change');

            case 6
                figure(iplot);
                for icfg=1:nconfig
                    meta=cfresults{icfg,1};
                    delstft=cfresults{icfg,2};
                    nframe=size(meta,1); % number of frames
                    plot((meta(1:nframe-1,1)'+meta(1:nframe-1,2)'-0.5)/fs,delstft(:,2));
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(delstft(:,2)))];
                    hold on
                end
                v_axisenlarge([-1 -1.05]);
                plot(get(gca,'xlim'),[0 0],':k');
                hold off
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('Inter-frame \Delta (dB)');
                legend(cfnamesmu,'location','best');
                title('Interframe Complex Change (Real-Imag interp)');

            case 7
                figure(iplot);
                for icfg=1:nconfig
                    meta=cfresults{icfg,1};
                    delstft=cfresults{icfg,2};
                    nframe=size(meta,1); % number of frames
                    plot((meta(1:nframe-1,1)'+meta(1:nframe-1,2)'-0.5)/fs,delstft(:,3));
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(delstft(:,3)))];
                    hold on
                end
                v_axisenlarge([-1 -1.05]);
                plot(get(gca,'xlim'),[0 0],':k');
                hold off
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('Inter-frame \Delta (dB)');
                legend(cfnamesmu,'location','best');
                title('Interframe Complex Change (Mag-Uphase interp)');

            case 8
                figure(iplot);
                for icfg=1:nconfig
                    meta=cfresults{icfg,1};
                    delstft=cfresults{icfg,2};
                    nframe=size(meta,1); % number of frames
                    plot((meta(1:nframe-1,1)'+meta(1:nframe-1,2)'-0.5)/fs,delstft(:,3)-delstft(:,2));
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(delstft(:,3)-delstft(:,2)))];
                    hold on
                end
                v_axisenlarge([-1 -1.05]);
                plot(get(gca,'xlim'),[0 0],':k');
                hold off
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('Inter-frame \Delta (dB)');
                legend(cfnamesmu,'location','best');
                title('Interframe Complex Change (Mag-Uphase minus Real-Imag)');

            case 9 %%% 8: single spectrum frequency (fplot) vs time
                figure(iplot);
                for icfg=1:nconfig
                    plot(cfresults{icfg,4}(:,1),db(abs(cfresults{icfg,4}(:,2))));
                    fminmax=cfresults{icfg,6}(2:3);
                    if fminmax(1)==fminmax(2)
                        cfnamesmu{icfg}=[cfnames{icfg} sprintf(' @ %.0f Hz',fminmax(1))];
                    else
                        cfnamesmu{icfg}=[cfnames{icfg} sprintf(' @ %.0f-%.0f Hz',fminmax)];
                    end
                    hold on
                end
                [spgbwt,spgbwf,spgbwb]=v_spgrambw(s,fs,'d'); % calculate spectrogram in dB
                [dum,fploti]=min(abs(spgbwf-fplot)); % find frequency index to plot
                dboff=20*log10(2*(length(spgbwt)-1)); % empirical gain correction
                plot(spgbwt,spgbwb(:,fploti)+dboff,':k');
                cfnamesmu{nconfig+1}=sprintf('Spectrogram @ %.0f Hz',spgbwf(fploti));
                v_axisenlarge([-1 -1.05]);
                plot(tplot*[1 1],get(gca,'ylim'),':k');
                hold off
                axlinkt=[axlinkt gca]; % link time axis
                xlabel('Time (s)');
                ylabel('STFT Power (dB)');
                legend(cfnamesmu,'location','best');
                title(sprintf('SFTF @ %.0f Hz',fplot));
                if figbold, figbolden; end;

            case 10 %%% 9: single spectrum frame (tplot) vs frequency
                figure(iplot);
                for icfg=1:nconfig
                    plot(cfresults{icfg,5}(:,1),db(abs(cfresults{icfg,5}(:,2))));
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' @ %.2f s',cfresults{icfg,6}(1))];
                    hold on
                end
                [spgbwt,spgbwf,spgbwb]=v_spgrambw(s,fs,'d'); % calculate spectrogram in dB
                [dum,tploti]=min(abs(spgbwt-tplot)); % find time index to plot
                dboff=20*log10(2*(length(spgbwt)-1)); % empirical gain correction
                plot(spgbwf,spgbwb(tploti,:)+dboff,':k');
                cfnamesmu{nconfig+1}=sprintf('Spectrogram @ %.2f s',spgbwt(tploti));
                v_axisenlarge([-1 -1.05]);
                plot(fplot*[1 1],get(gca,'ylim'),':k');
                hold off
                axlinkf=[axlinkf gca]; % link time axis
                xlabel('Freq (Hz)');
                ylabel('STFT Power (dB)');
                legend(cfnamesmu,'location','best');
                title(sprintf('SFTF @ %.2f s',tplot));
                if figbold, figbolden; end;

            case 11 %%% 10: mean changes by phone type (complex: mag-phase interp)
                figure(iplot);
                for icfg=1:nconfig
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(cfresults{icfg,2}(:,3)))];
                    hold on
                end
                bar(delstftmu(:,:,3));
                set(gca,'xtick',1:nphntype,'xticklabel',phntypename);
                ylim=get(gca,'ylim');
                legend(cfnamesmu,'location','best');
                ylabel('\Delta (dB)')
                title('Interframe Complex Change (Mag-Uphase interp)');
                v_fig2pdf(gcf,figpdf);

            case 12 %%% 11: mean changes by phone type (magnitude)
                figure(iplot);
                for icfg=1:nconfig
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(cfresults{icfg,2}(:,1)))];
                    hold on
                end
                bar(delstftmu(:,:,1));
                set(gca,'xtick',1:nphntype,'xticklabel',phntypename);
                legend(cfnamesmu,'location','northeast');
                ylabel('\Delta (dB)')
                title('Interframe Magnitude Change');
                if dopdf, v_fig2pdf(figpdf), end;

            case 13 %%% 12: mean changes by phone type (complex: real-imag interp)
                figure(iplot);
                for icfg=1:nconfig
                    cfnamesmu{icfg}=[cfnames{icfg} sprintf(' \\mu=%.1f',mean(cfresults{icfg,2}(:,2)))];
                    hold on
                end
                bar(delstftmu(:,:,2));
                set(gca,'xtick',1:nphntype,'xticklabel',phntypename);
                legend(cfnamesmu,'location','best');
                ylabel('\Delta (dB)')
                title('Interframe Complex Change (Real-Imag interp)');
                if dopdf, v_fig2pdf(figpdf), end;

            case 14 %%% 13: within-phone Magnitude Variance
                figure(iplot);
                for icfg=1:nconfig
                    cfnamesmu{icfg}=cfnames{icfg};
                    hold on
                end
                bar(stftmpdbmu(:,:,1));
                set(gca,'xtick',1:nphntype,'xticklabel',phntypename);
                legend(cfnamesmu,'location','northeast');
                ylabel('\Delta (dB)')
                title('within-phone Magnitude Variance');
                if dopdf, v_fig2pdf(figpdf), end;

            case 15 %%% 14: within-phone Complex Variance
                figure(iplot);
                for icfg=1:nconfig
                    cfnamesmu{icfg}=cfnames{icfg};
                    hold on
                end
                bar(stftmpdbmu(:,:,2));
                set(gca,'xtick',1:nphntype,'xticklabel',phntypename);
                legend(cfnamesmu,'location','best');
                ylabel('\Delta (dB)')
                title('within-phone Complex Variance');
                if dopdf, v_fig2pdf(figpdf), end;
        end
    end
end
v_tilefigs;
if ~isempty(axlinkt)
    linkaxes(axlinkt,'x');
end
if ~isempty(axlinkf)
    linkaxes(axlinkf,'x');
end
if ~isempty(axlinknm)
    linkaxes(axlinknm,'x');
end
