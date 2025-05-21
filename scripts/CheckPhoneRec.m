% TIMIT phone recognition with Gaussian mixtures
%
% Ideas:
% (1) Try varying the number of mel coefficients
% (2) Train and/or test only on the phone centres (i.e. ignore transition regions)
% (3) Use LDA on log mel spectrum instead of DCT
% (4) Incorporate transition probabilities
clear;
close all;
%
parvar={'cfname' 'epoch' 'nmix'};                                                              % cell array row listing the parameters to change during trials
configpars={
{'30ms-g1' 0 1};   % one row per trial giving the parameter values
    % {'30ms-m4' 0 4};
    {'ep-g1' 1 1};
    {'ep-g1' 1 1};
   {'epq-natm-g1' 1 1 'smoothalg' 'quadlog' 'interpstft' 'natural' 'interpdom' 'magcph'};
    % {'ep-g4' 1 4};
    {'ep-g1-lda15' 1 1 'ldatype' 'lda'};
    % {'e-g4-lda15'  1 4 'ldatype' 'lda'};
    };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plots (include in plotlist to plot, * = final test utterance only):
%           1       List of configurations
%           7       Accuracy bar graph
%     *   500       Spectrogram of final test utterance
%    ... plots for each configuration (added to configuration number).
%         100+      Between/Within variances for each feature element
%         200+      Mean and stddev of each feature for each phone
%         300+      Recognition performance by phone
%         400+      Phone confusion matrix
%     *   600+      Feature-Spectrogram + recognition results
%     *   800+      Mel-Spectrogram + recognition results
%
parplt=[1:10]; 								% configurations to plot
pltsel=[1 7 300:100:800]; 				    % list of plots to do
nfile=[100 100];                      		% number of TIMIT files for training and testing
% nfile=[10 10];                      		% number of TIMIT files for training and testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Default parameter values     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% = global =
par0.toptail=2;                      % 0 = keep all frames, 1 = Remove leading/training silence, 2 = keep only vowels
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
% = stftgrid =
par0.interpstft=  'none';            % interpolation method for call to griddata: {'none','nearest','linear','natural','cubic','v4'}
par0.interpgrid= [];                 % par.interpgrid=[nhop nbin] is calculated below from par.interpflen, par.interpov and fs
par0.interpflen=  0.007;             % target frame length in seconds (reciprocal of frequency-grid spacing) [0.01]
par0.interpov=    1;                 % overlap factor; 1 for no overlap, 2 for 50%; time grid-spacing is par.interpflen/par.interpov
par0.interpbph=   0.196;             % time-frequency distance tradeoff in bins/hop
par0.interpdom=   'cplx';            % interpolation domain: {'cplx','magcph','crmcph'}
% = stfte2melc
par0.nmel=        29;                % number of mel bins = 29 @ 16kHz
par0.fbank=       'm';               % filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel}
par0.keepDC=      0;                 % preserve DC as the lowest MEL bin {0, 1}
par0.MELphase=    'linear';          % MEL STFT phase calculation: {'true','zero','linear','piecewiselin','piecewiselind','grpdel'}
par0.MELdom=      'mag';             % MEL filterbank domain: {'mag', 'pow'}
par0.regwt=       0.01;              % regularizing factor for phase weights when par.MELphase='piecewiselin'
par0.loops=       5;                 % maximum number of iterations when par.MELphase='piecewiselin'
par0.ssqthr=      1.0;               % stopping threshold when par.MELphase='piecewiselin'
% = melc2stfte
par0.invmethod=  'pinv';             %  MEL inversion method: {'pinv','transp','linterp'} ['linterp']
% = mel2melcep
par0.ncep=12;                                % number of cepstral coefficients
par0.wcep='d0';                              % flags for delta and delta-delta coefficients
% = LDA feature compression
par0.nlda=15;                                % number of parameters after LDA (or 0 for keep same dimension)
par0.ldatype='none';                         % LDA algorithm: {'none','lda','hda'}
% = phone recognition
par0.ndp=0;                                  % number of previous paths to keep in Dynamic programming (0 for no DP,>=67 for all phones)
par0.tranwt=0.6;                             % weighting for transition log probability (0=ignore transition probs, 1=take full account)
par0.fullcov=0;                              % Covariance matrices: 0=diagonal, 1=full
par0.nmix=1;                                 % number of mixtures per phone
%
% set global parameters
%
fs=16000;                                   % TIMIT sample frequency
timflg={'n','t'};                           % TIMIT file options for training and testing
%
% derived parameters
%
ncfg=size(configpars,1);                        % number of parameter configurations
% nmel=floor(3*log(fs));                      % number of mel bins = 29 @ 16kHz
[sph,vph]=w_phoncode('tU');                 % list of all TIMIT phones
nphons=length(sph);                         % number of possible TIMIT phones [67]
silphon=64;                                 % silence phone number (currently 64)
ty=cell(1,2);                               % initialize list of TIMIT files to use
for itt=1:2 % training and testing
    ty{itt}=timitfiles(timflg{itt},nfile(itt)); % list of TIMIT files to use
end
% initialize per-conficuration statistics
cfnames=cell(ncfg,1);                                    % list of configuration names
cfresults=cell(ncfg,1);                                  % space for meta information
acc=zeros(1,ncfg); % overall accuracy
for icfg=1:ncfg                                         % loop through parameter configurations
    [par,partxt]=parupdate(par0,configpars{icfg},parvar,'%cfname%: $c&=%&%$, $'); % Set up the parameters for this trial
    cfnames{icfg}=par.cfname;                                   % save configuration name
    %
    % derived parameters that depend on parameter configuration
    %
    needlda=~strcmp(par.ldatype,'none');
    ndp=par.ndp; % number of alternatives to keep in the dynamic programming
    nd=(par.ncep+any(par.wcep=='0'))*(1+any(par.wcep=='d')+any(par.wcep=='D'));    % size of raw feature vector
    np=nd-(needlda && par.nlda>0 && par.nlda<nd)*(nd-par.nlda);                    % actual feature dimension: either nd or par.nlda
    needfullcov=par.fullcov || needlda; % need to calculate full covariance matrices from the training data
    %
    if needfullcov
        ndr=repmat(1:nd,1,nd);                          % row indices for nd*nd covariance matrix
        ndc=floor((0:nd^2-1)/nd)+1;                     % col indices for nd*nd covariance matrix
        ndd=1:nd+1:nd^2;                                % diagonal element indices
    else
        ndr=1:nd;
        ndc=1:nd;
        ndd=1:nd;                                       % diagonal element indices
    end
    nfft=round(par.fixfl*fs);                           % frame length in samples [~ par.fixfl seconds]
    inc=floor(nfft/2);                                  % frame increment = 50% of frame length
    for itt=1:2                                         % training and testing
        if itt==1                                       % training initialization
            vn=zeros(nphons,1);                         % initialize: phone count
            vm=zeros(nphons,nd);                        % ... phone sum
            vc=zeros(nphons,length(ndc));               % ... phone sum of squares
            featall=[];                                 % ... feature list
            phvkxtall=[];                               % ... phone label list
            transm=sparse(nphons,nphons);               % ... transition counts
        else                                            % testing initialization (done at the start of pass 2)
            vv=vc(:,ndd);                               % extract diagonal elements from full or diagonal sums
            % vn(silphon)=0;                            % eliminate the silence phone from training data
            % vm(silphon,:)=0;
            % vv(silphon,:)=0;
            mk=vn>0;                                    % mask for phones that occur in training data
            phonlist=find(mk);                          % list of valid phones that had training data
            phonlistsil=find(phonlist==silphon);        % silence phone in list (should always exist)
            vnk=vn(mk);                                 % occurrence frequencies of valid phones (used as GMM weights when normalized)
            gvn=length(vnk);                            % number of valid phone classes
            vph=sum(vnk);                               % total number of valid training frames
            vnkw=vnk/vph;                               % phone weights (sum to unity)
            vphm=sum(vm,1)/vph;                         % global mean feature vector
            vphv=sum(vv,1)/vph-vphm.^2;                 % global feature vector variance
            vphs=sqrt(vphv);                            % global feature vector std dev
            vmk=vm(mk,:)./repmat(vnk,1,nd);             % calculate phone means
            vvk=vv(mk,:)./repmat(vnk,1,nd)-vmk.^2;      % calculate within-phone variances
            if needfullcov
                vck=vc(mk,:)./repmat(vnk,1,nd^2)-vmk(:,ndr).*vmk(:,ndc);        % calculate within-phone covariances
                vck(:,1:nd+1:nd^2)=vck(:,1:nd+1:nd^2)+repmat(vphv*1e-6,gvn,1);  % add fraction of global variance to ensure positive definiteness
                vck=reshape(vck',nd,nd,gvn);            % reshape so phone number is the last dimension
            else
                vck=vvk;                                % within-phone diagonal covariances
            end
            gvm=mean(vmk,1);                            % unweighted global mean i.e. mean of phone means
            gvv=sum(vmk.^2,1)/gvn-gvm.^2;               % global between-phone variance
            gvc=vmk'*vmk/gvn-gvm'*gvm;                  % between-phone covariance matrix
            switch par.ldatype
                case 'none'
                    ldam=diag(1./vphs); % variance-normalizaing transformation
                case 'lda'
                    wcc=mean(vck,3);                        % unweighted mean within-phone covariance matrix
                    ldam=v_dualdiag(wcc,gvc);
                    ldam=ldam(:,1:np);                      % only keep the first np columns
                case 'hda'
                    wcc=mean(vck,3);                        % unweighted mean within-phone covariance matrix
                    ldam=v_dualdiag(wcc,gvc);               % do standard LDA first
                    ldam=ldam(:,1:np);                      % only keep the first np columns
                    opt = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','notify');
                    f=@(a)hdadobj(a,vmk,vck,vnkw);
                    ldam=fminunc(f,ldam,opt);
            end
            vmpk=vmk*ldam;
            if par.fullcov                                  % par.fullcov=1
                vcpk=zeros(np,np,gvn);
                for iph=1:gvn
                    vcpk(:,:,iph)=ldam'*vck(:,:,iph)*ldam;
                end
            elseif needfullcov                              % par.fullcov=0 && needlda
                vcpk=zeros(gvn,np);
                for iph=1:gvn
                    vcpk(iph,:)=diag(ldam'*vck(:,:,iph)*ldam)';
                end
            else                                            % par.fullcov=0 && ~needlda
                vcpk=vck.*repmat(diag(ldam)'.^2,gvn,1);
            end
            % calculate the GMMs
            featall=featall*ldam;                       % transform the data using the LDA matrix
            nmix=par.nmix;                              % number of gmm mixtures
            gmmm=zeros(nmix,np,gvn); % GMM means
            gmmw=zeros(nmix,gvn); % GMM weights
            if par.fullcov                                  % par.fullcov=1
                gmmv=zeros(np,np,nmix,gvn);  % GMM full covariance matrices
                for iph=1:gvn
                    pmsk=(phvkxtall==phonlist(iph)); % mask for frames corresponding to this phoneme
                    [gmmm(:,:,iph),gmmv(:,:,:,iph),gmmw(:,iph)]=v_gaussmix(featall(pmsk,:),[],[],nmix,'v');    % create GMM with k mixtures and full covariances
                end
                % [diag(gmmv(:,:,1,1))'; diag(vcpk(:,:,1))'; diag(gmmv(:,:,1,1))'-diag(vcpk(:,:,1))'] % compare phone #1 ********** DEBUG ONLY ***************
            else
                gmmv=zeros(nmix,np,gvn);  % GMM diagonal covariance matrices
                for iph=1:gvn
                    pmsk=(phvkxtall==phonlist(iph)); % mask for frames corresponding to this phoneme
                    [gmmm(:,:,iph),gmmv(:,:,iph),gmmw(:,iph)]=v_gaussmix(featall(pmsk,:),[],[],nmix);    % create GMM with k mixtures and full covariances
                end
                % [gmmv(1,:,1); vcpk(1,:); gmmv(1,:,1)-vcpk(1,:)] % compare phone #1 ********** DEBUG ONLY ***************
            end
            % Calculate F ratio as the trace of (W\B) where B and W are the diagonal between and average within-phone variances
            % should maybe weight within-phone variances by the frequency of the phones but I'm not sure that others do this.
            wcv=mean(vvk,1);                            % mean within-phone variance
            frat=sum(gvv./wcv);                         % generalized F ratio
            confm=sparse(nphons,nphons);                % initialize confusion matrix
            if ndp>0
                transmk=1+transm(mk,mk);
                translp=par.tranwt*log(transmk./repmat(sum(transmk,2),1,gvn)); % log transition probs (add 1 to pevent -Inf values)
            end
        end
        nfiles=length(ty{itt});
        for ifile=1:nfiles                                          % loop through training or testing files
            v_finishat([icfg 1 ncfg; itt 1 2; ifile 1 nfiles]);     % estimate finish time
            [s,fs,wrd,phn]=gettimit(ty{itt}{ifile},'n');            % read timit file and normalize power
            [uniph,phv]=w_phoncode('tFU',phn(:,2),11);              % convert phone symbols to TIMIT phone numbers
            ophnt=cell2mat(phn(:,1));                               % start and end times of phones [seconds]
            nphone=size(ophnt,1);                                   % number of phones in this utterance
            if par.epoch                                            % *** epoch-based frames ***
                frameend=gs_frames(s,fs,par);                       % last sample of each frame
                framelim=[[1 frameend(1:end-1)+1];frameend];        % create frame limits
                framelen=1+[-1 1]*framelim;                         % length of each frame in samples
                metain=[framelim(1,:);framelen]';
                [stft,meta]=stfte(s,metain,[],par);                 % epoch-based STFT
                [stft,meta]=stftgrid(stft,meta,par);                % optionally map onto a fixed grid
                [melbm,melbu,cfhz,mbm]=stfte2melc(stft,fs,par.nmel,par);
                c=mel2melcep(melbm,par.wcep,par.ncep);
                tc=(meta(:,1:2)*[1;0.5]-1.5)/fs;                    % frame centres in seconds (1st sample @ zero)
                pth=max(melbm(:))*1E-20;                            % low threshold to avoid infinite logs
                y=log(max(melbm,pth));                              % log mel spectrogram
            else                                                    % *** fixed frames ***
                [z,tc]=v_enframe(s,0.54-0.46*cos(2*pi*(0:nfft-1)'/(nfft-1)),inc); % Hamming window
                tc=(tc-1)/fs;                                       % times of frame centres with 1st sample @ 0 [seconds]
                ss=v_rfft(z,nfft,2);
                [c,tcmc,y,mc]=w_melcepstx(ss,fs,[par.wcep 'S'],par.ncep,par.nmel,nfft);       % do mel cepstrum with input as stft
            end
            nframe=size(c,1);                                       % number of frames
            [tt,ixt,jxt]=v_sort([tc;ophnt(:,1)]); % log mel spectrum
            kxt=min(max(jxt(1:nframe)-(1:nframe)',1),nphone);       % phone number of each frame
            if par.toptail                                          % check if need to remove initial and final silence segment
                switch par.toptail
                    case 1                                          % par.toptail=1: remove leading/trailing silent frames
                        phvkxt=phv(kxt,1);                          % TIMIT phone number for each frame
                        toptail=[1:find(phvkxt~=silphon,1)-2 find(phvkxt~=silphon,1,'last')+2:nframe]; % leading and trailing silent frames
                    case 2
                        toptail=phv(kxt,3)~=14;                     % keep only vowel frames
                end
                c(toptail,:)=[];                                    % remove unwanted frames from cepstrum
                y(toptail,:)=[];                                    % ... and log mel spectrum
                tc(toptail)=[];                                     % ... and frame centre times
                kxt(toptail)=[];                                    % ... and original phone sequence number of each frame
                nframe=size(c,1);                                   % update number of frames
            end
            phvkxt=phv(kxt,1);                                      % TIMIT phone number for each frame
            if itt==1                                               % *** training ***
                fmap=sparse(phvkxt,1:nframe,1,nphons,nframe);       % matrix to map frames to TIMIT phones
                transm=transm+sparse(phvkxt(1:end-1),phvkxt(2:end),1,nphons,nphons); % accumulate transition counts
                featall=[featall;c];                                % add features onto the global list
                phvkxtall=[phvkxtall; phvkxt];                      % add phone labels onto the global list
                vn=vn+sum(fmap,2);                                  % accumulate phone counts
                vm=vm+fmap*c;                                       % accumulate coefficient sums
                vc=vc+fmap*(c(:,ndr).*c(:,ndc));                    % accumulate full covariances
            else         % *** testing ***
                % [lp,rp,kh]=v_gaussmixp(c*ldam,vmpk,vcpk,vnkw);         % phonlist classification for each frame
                lp=zeros(nframe,gvn); % space for classifier log probabilities
                if par.fullcov
                    for iph=1:gvn
                        lp(:,iph)=v_gaussmixp(c*ldam,gmmm(:,:,iph),gmmv(:,:,:,iph),gmmw(:,iph));         % phonlist classification for each frame
                    end
                else
                    for iph=1:gvn
                        lp(:,iph)=v_gaussmixp(c*ldam,gmmm(:,:,iph),gmmv(:,:,iph),gmmw(:,iph));         % phonlist classification for each frame
                    end
                end
                if ndp>0
                    if ndp<gvn
                        xph=repmat(phonlistsil,ndp,1);              % previous phone (initialize to silence)
                        xutil=zeros(ndp,1);                         % initialize DP costs
                        xprev=ones(nframe,ndp);                     % traceback
                        xphon=ones(nframe,ndp);                     % best ndp phone choices for each frame
                        for ifr=1:nframe
                            yutil=repmat(xutil,1,gvn)+translp(xph,:)+repmat(lp(ifr,:),ndp,1); % cumulative costs
                            [phutil,iprv]=max(yutil,[],1);          % find the best choice for each value of current phone
                            [zutil,iy]=sort(phutil,'descend');      % Find the best utilities of current frame for each phone
                            xutil=zutil(1:ndp)';                    % ... but only keep the best few utilities
                            xph=iy(1:ndp)';                         % ... and corresponding phone choice
                            xphon(ifr,:)=xph;                       % best ndp phone choices for the current frame
                            xprev(ifr,:)=iprv(xph);                 % ... and corresponding traceback indices
                        end
                        % now do traceback
                        ilst=1;                                     % list choice for the final frame is always 1
                        kh(nframe)=iy(ilst);                        % final frame is the phone with the highest utility
                        for ifr=nframe-1:-1:1
                            ilst=xprev(ifr+1,ilst);                 % list entry for traceback
                            kh(ifr)=xphon(ifr,ilst);                % phon for this frame
                        end
                    else                                            % keep all phones as possibilities
                        xph=repmat(phonlistsil,gvn,1);              % previous phone (initialize to silence)
                        xutil=zeros(gvn,1);                         % initialize DP costs
                        xprev=ones(nframe,gvn);                     % traceback
                        xphon=ones(nframe,gvn);                     % best gvn phone choices for each frame
                        for ifr=1:nframe
                            yutil=repmat(xutil,1,gvn)+translp(xph,:)+repmat(lp(ifr,:),gvn,1); % cumulative costs
                            [phutil,iprv]=max(yutil,[],1);          % find the best choice for each value of current phone
                            xutil=phutil';                          % becomes the previous costs for next frame
                            xprev(ifr,:)=iprv;                      % ... and corresponding phone for previous frame
                        end
                        % now do traceback
                        [tutil,ilst]=max(xutil);                    % maximize utility for final frame
                        kh(nframe)=ilst;                            % final frame is the phone with the highest utility
                        for ifr=nframe-1:-1:1
                            ilst=xprev(ifr+1,ilst);                 % list entry for traceback
                            kh(ifr)=ilst;                           % phon for this frame
                        end
                    end
                else
                    [lpdum,kh]=max(lp,[],2);                        % if no DP, just take maximum log prob for each frame
                end
                confm=confm+sparse(phvkxt,phonlist(kh),ones(nframe,1),nphons,nphons); % confusion matrix
            end
        end                                                         % loop through files
    end                                                             % training/testing
    % ******** calculate statistics
    testcnt=full(sum(confm,2));                                     % number of test phones in each phone-class
    corrcnt=full(diag(confm));                                      % number of correctly recognised phones of each phone-class
    normcf=sparse(1:gvn,1:gvn,testcnt(mk).^-1);                     % diagonal count-normalizing matrix
    acc(icfg)=sum(corrcnt(mk))/sum(testcnt(mk));                    % accuracy: fraction of phones correctly recognised
    % ******** now do configuration-dependent plotting
    if any(icfg==parplt)                                            % need to do plots for this configuration
        if any(pltsel==100)                                           % 100: Between/Within variances for each feature element
            figure(100+icfg);
            bar(gvv./wcv);                              % between-phone to within-phone variance ratio for each feature
            v_texthvc(0.95,0.95,sprintf('%d Train Files, %d Frames\n%d Phones\nRatio Sum = %.2f',length(ty{1}),vph,gvn,frat),'RTk');
            xlabel('Feature Element');
            ylabel('Between / Within Variance');
            title(sprintf('Training: %s',par.cfname));
        end
        if any(pltsel==200)                               % 200: Mean and stddev of each feature for each phone
            figure(200+icfg);
            subplot(2,1,2);
            if par.fullcov
                vcpkd=reshape(vcpk,nd^2,gvn); % reshape to have vectorized covariance in each column
                imagesc(sqrt(vcpkd(ndd,:))); % extract diagonal elements
            else
                imagesc(sqrt(vcpk)');
            end
            colorbar;
            axis('xy');
            if ~needlda
                hold on; plot([0.5 gvn+0.5],[12.5 12.5],'w-'); hold off;
            end
            set(gca,'xtick',1:gvn,'xticklabel',sph(mk));
            xlabel('Phone');
            ylabel('Std Dev');
            subplot(2,1,1);
            imagesc(vmpk'); % normalize the means
            colorbar;
            axis('xy');
            if ~needlda
                hold on; plot([0.5 gvn+0.5],[12.5 12.5],'w-'); hold off;
            end
            set(gca,'xtick',1:gvn,'xticklabel',sph(mk));
            ylabel('Mean');
            title(sprintf('Training: %s',par.cfname));
        end
        if any(pltsel==300) %  300: Recognition performance by phone
            figure(300+icfg);
            % bar([corrcnt(mk) testcnt(mk)-corrcnt(mk)],'stacked');
            % legend('Correct','Error','location','northwest');
            bar(100*corrcnt(mk)./(testcnt(mk)+(testcnt(mk)==0)));
            hold on
            plot([0.5 gvn+0.5],[100 100]*acc(icfg),'--k');
            hold off
            set(gca,'xtick',1:gvn,'xticklabel',sph(mk),'ylim',[0 100]);
            v_texthvc(0.98,0.98,sprintf('%d Test Files, %d Frames\n%.0f%% Correct',length(ty{2}),sum(testcnt(mk)),100*acc(icfg)),'RTk');
            % set(gca,'xtick',1:gvn,'xticklabel',sph(mk),'yscale','log','ylim',max(testcnt(mk)).^[-0.05 1.05]);
            % yticklab=get(gca,'yticklabel');
            % simult={'' 'k' 'M' 'G'};
            % for i=1:length(yticklab)
            %     ylabi=yticklab{i};
            %     if length(ylabi)>3 && strcmp(ylabi(1:4),'10^{')
            %         ylexp=str2num(ylabi(5:end-1)); % find exponent
            %         ylexp3=min(max(floor(ylexp/3),0),length(simult)-1);
            %         yticklab{i}=[num2str(10^(ylexp-3*ylexp3)) simult{ylexp3+1}];
            %     end
            % end
            % set(gca,'yticklabel',yticklab);
            xlabel('Phone');
            ylabel('% correct');
            title(sprintf('Reconition: %s',par.cfname));
        end
        if any(pltsel==400) %  400: Confusion Matrix
            figure(400+icfg);
            imagesc(confm(mk,mk)'*normcf); %%%%%% we should include all phones in test *or* training sets
            hold on;
            dol=(1:2*gvn+1)'/2;
            plot([floor(dol)+0.5 ceil(dol)-0.5],[ceil(dol)-0.5 floor(dol)+0.5],'-w');
            hold off;
            set(gca,'xtick',1:gvn,'xticklabel',sph(mk));
            set(gca,'ytick',1:gvn,'yticklabel',sph(mk));
            ylabel('Estimated');
            xlabel('True');
            title(sprintf('Confusion: %s',par.cfname));
        end
        if any(pltsel==600) %  600: Cepstrogram + recognition results
            figure(600+icfg);
            % phcor=(2*(phv(kxt,1)==phonlist(kh))-1).*(phv(kxt,1)~=silphon); % frame recognition: -1=wrong, 0=silence, 1=right
            phcor=(2*(phv(kxt,1)==phonlist(kh))-1); % frame recognition: -1=wrong, 1=right
            phcs=sum([phcor abs(phcor)]);
            phacc=50*sum(phcs)/phcs(2); % frame accuracy
            % cnorm=(c-repmat(vphm,nframe,1))./repmat(vphs,nframe,1);
            cnorm=c*ldam;
            imagesc(1:nframe,0:np,[max(abs(cnorm(:)))*phcor cnorm]');
            axis xy;
            if ~needlda
                hold on;
                plot([1; nframe],[1; 1]*[0.5 (par.ncep+0.5)],'-w');
                hold off;
            end
            title(sprintf('%s: Normalized Features',par.cfname));
            ylabel('Feature');
            ax1=gca;
            set(ax1,'XaxisLocation','top');
            kxtd=find(kxt(1:end-1)~=kxt(2:end)); % last frame of each true phone (except the last one)
            hold on
            plot([1;1]*(kxtd'+0.5),get(ax1,'ylim')','-k')
            hold off
            ytickl=get(ax1,'yticklabel');
            set(ax1,'ytick',get(ax1,'ytick'),'yticklabel',[{sprintf('%.0f%%',phacc)}; ytickl(2:end)]); % relace '0' label with % score
            ax2=axes('Position',ax1.Position,'XaxisLocation','bottom','YAxisLocation','right','Color','none','YTick',[],'TickLength',[.00014 .0014]);
            set(ax2,'xlim',[0.5 nframe+0.5],'Xtick',[[1;kxtd+1] [kxtd;nframe]]*[0.5;0.5],'Xticklabel',uniph(kxt([kxtd;nframe])));
            % [tdum,idum,jjx]=v_sort([ophnt*[0.5;0.5];tc]);
            % jx=(jjx(1:nphone)-(1:nphone)')/nframe; % phone centres
            % jxmsk=jx>0 & jx<1; % mask for phones that lie within the axis limits
            % set(ax2,'Xtick',jx(jxmsk),'Xticklabel',uniph(jxmsk));
            xlabel('Frame #, True Phone');
        end
        if any(pltsel==800) %  800: Mel-Spectrogram + recognition results
            figure(800+icfg);
            % phcor=(2*(phv(kxt,1)==phonlist(kh))-1).*(phv(kxt,1)~=silphon); % frame recognition: -1=wrong, 0=silence, 1=right
            phcor=(2*(phv(kxt,1)==phonlist(kh))-1); % frame recognition: -1=wrong, 1=right
            phcs=sum([phcor abs(phcor)]);
            phacc=50*sum(phcs)/phcs(2); % frame accuracy
            imagesc(1:nframe,0:size(y,2)-1,[max(abs(y(:)))*phcor y]');
            axis xy;
            title(sprintf('%s: Logmel Spgram',par.cfname));
            ylabel('Mel bin');
            ax1=gca;
            set(ax1,'XaxisLocation','top');
            kxtd=find(kxt(1:end-1)~=kxt(2:end)); % last frame of each true phone (except the last one)
            hold on
            plot([1;1]*(kxtd'+0.5),get(ax1,'ylim')','-k')
            hold off
            ytickl=get(ax1,'yticklabel');
            set(ax1,'ytick',get(ax1,'ytick'),'yticklabel',[{sprintf('%.0f%%',phacc)}; ytickl(2:end)]); % relace '0' label with % score
            ax2=axes('Position',ax1.Position,'XaxisLocation','bottom','YAxisLocation','right','Color','none','YTick',[],'TickLength',[.00014 .0014]);
            set(ax2,'xlim',[0.5 nframe+0.5],'Xtick',[[1;kxtd+1] [kxtd;nframe]]*[0.5;0.5],'Xticklabel',uniph(kxt([kxtd;nframe])));
            xlabel('Frame #, True Phone');
        end
    end
    cfresults{icfg,1}=partxt;                            % text description of configuration
end % loop through parameter configurations
%
% configuration-independent or summary plota
%
if any(pltsel==1)                                       %  1: List of configurations
    figure(1);
    cfdesc='';
    linelen=40;
    for i=1:ncfg
        cfdesci=sprintf('%d) %s',i,cfresults{i,1});
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
end
if any(pltsel==500)                                     %  500: Spectrogram
    figure(500);
    phnu=phn;
    phnu(:,2)=cellstr(uniph);
    v_spgrambw(s,fs,'pJcwatmf',[],[],[],[],{wrd phnu});
    title('Test Sentence');
end
if any(pltsel==7)                                       %  7: Accuracy
    figure(7);
    bar(100*acc);
    set(gca,'xtick',1:ncfg,'xticklabel',cfnames);
    ylabel('% Accuracy');
    xlabel('Configuration');
end
v_tilefigs;