% test copy synthesis:  signal -> stfte -> melc -> stfte -> signal
close all;
% Define default parameters
nmel=30;
wavopt=2;
viamel=0; % 0= stfte/istfte only, 1=stfte/stfte2melc/melc2stfte/istfte
%
par.pitchlim=[40 50 400];                   % {min targt max} pitch (Hz)
par.gcifrac=0.3;                            % position of GCI in analysis frame
par.GCImethod='YAGA';                       % Glottal closure detection algorithm
par.offset='none';                          % offset removal: {'none','mean'}
par.scale='none';                           % scaling method: {'none','peakabs','rms'}
par.pad='zero';                             % zero-padding method: {'none','zero','ends'}
par.groupdelay='none';                      % linear phase component: {'none','ewgd','ewgdint','cplx','cplxint'}
par.keepDC=1;                               % preserve DC as the lowest MEL bin {0, 1}
par.MELphase='piecewiselin';                % MEL STFT phase calculation: {'true','zero','linear','piecewiselin'}
par.MELdom='pow';                           % MEL filterbank domain: {'mag', 'pow'}
par.keepDC=1;                               % preserve DC in the lowest MEL bin {0 or 1}
par.invmethod='linterp';                    % MEL inversion method: {'pinv','transp','linterp'}
par.MELphase='piecewiselin';                % MEL STFT phase calculation: {'true','zero','linear','piecewiselin'}
par.MELdom='pow';                           % MEL filterbank domain: {'mag', 'pow'}
par.regwt=0.01;                             % regularizing factor for phase weights when par.MELphase='piecewiselin'
par.loops=5;                                % maximum number of iterations when par.MELphase='piecewiselin'
par.ssqthr=1.0;                             % stopping threshold when par.MELphase='piecewiselin'
%
%
% create the data waveform: must set s, fs, framelim, frameend
%
switch wavopt
    case 0 % keep previous
    case 1 % cosine wave
        nper=10; % period of cosine wave
        frameend=round(cumsum([20 20 30 40 50 60 20.14 30.14 40.14])*nper);
        framelim=[[1 frameend(1:end-1)+1];frameend];  % create frame limits
        framelen=1+[-1 1]*framelim;                 % length of each frame in samples
        % create cosine wave
        smax=2;
        s=smax*cos(2*pi*(1:framelim(end))'/nper);
        fs=8000;
    case 2 % glottal flow derivative
        fs=16000; % sample frequency
        periods=[10 10 8 6 4 6 8 10]*1e-3; % periods in ms
        frac0=0; % fraction of a sample at each frame boundary
        framelen=round(periods*fs);
        frameend=cumsum(framelen);
        framelim=[[1 frameend(1:end-1)+1];frameend];  % create frame limits
        nframe=length(framelen);
        ns=frameend(end);
        ts=zeros(ns,1); % space for times
        for i=1:nframe
            ts(framelim(1,i):framelim(2,i))=i-1+frac0+(0:framelen(i)-1)/framelen(i);
        end
        s=v_glotlf(1,ts);
    case 3
        [timf,ttimy,tk,s,fs,timwrd,timphn]=timitfiles('pz',1);    % read a level-normalized random TIMIT file
        frameend=gs_frames(s,fs,par);
        framelim=[[1 frameend(1:end-1)+1];frameend];  % create frame limits
        framelen=1+[-1 1]*framelim;                 % length of each frame in samples
end
parvar={'pad' 'invmethod'};                 % cell array row listing the parameters to change during trials
partry={'zero' 'linterp'}; % one row per trial giving the parameter values
metain=[framelim(1,:);framelen]';
for i=1:size(partry,1)
    %
    %%%%%%%%%%% First set up the parameters for this trial %%%%%%%%%%%
    %
    partxt='';
    for j=1:length(parvar)
        par.(parvar{j})=partry{i,j};
        if ischar(partry{i,j})
            partxt=[partxt ', ' parvar{j} '=' partry{i,j}];
        else
            partxt=[partxt ', ' parvar{j} '=' num2str(partry{i,j})];
        end
    end
    partxt(1:2)=[];
    %
    %%%%%%%%%%%%% Now do the processing %%%%%%%%%%%%%%%%%%
    %
    [stft,meta]=stfte(s,metain,[],par);
    if viamel
        [melbm,melbu,cfhz,mbm]=stfte2melc(stft,fs,nmel,par);
        stftr=melc2stfte(melbm,melbu,fs,meta,par);
    else
        stftr=stft; % copy stft directly
    end
    sr=istfte(stftr,meta);
    %
    %%%%%%%%%%%%%%%% now do plots %%%%%%%%%%%%%%%%%%
    %
    figure(10+i);                                           % ********* plot input signal
    plot(s);
    hold on
    fredge=unique([meta(:,1)-0.5; meta(:,1:2)*[1;1]-0.5]);  % frame boundaries
    if ~all(meta(:,4)==0)                                   % if there are any non-zero offsets
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[1;1]*meta(:,4)',':k');
    end
    if ~all(meta(:,5)==1)                                   % if there are any non-unity scale factors
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[1;1]*meta(:,5)','--r');
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[-1;-1]*meta(:,5)','--r');
    end
    v_axisenlarge([-1 -1.05]);
    ylim=get(gca,'ylim');
    plot([1;1]*fredge',ylim','k:');                         % mark frame boundaries
    if ~all(meta(:,6)==0)                                   % if there are any non-zero group delays
        plot(meta(:,[1 6])*[1;1]-1,ylim*[0.97;0.03],'^r',meta(:,[1 6])*[1;1]-1,ylim*[0.03;0.97],'vr');
    end
    hold off
    xlabel('Sample Number');
    ylabel('Input Signal');
    title(['Input: ' partxt]);
    figure(20+i);                                           % ********* plot stfte
    stfte(s,metain,[],par);
    title(['STFTE in: ' partxt]);
    if viamel
        figure(30+i);                                           % ********* plot mel spectrum
        stfte2melc(stft,fs,nmel,par);
        title(['Mel: ' partxt]);
        figure(40+i);                                           % ********* plot output stfte
        melc2stfte(melbm,melbu,fs,meta,par);
        title(['STFTE out: ' partxt]);
    end
    figure(60+i);                                           % ********* plot output signal
    plot(sr,'-b');
    hold on
    plot(sr-s,'-r');
    v_axisenlarge([-1 -1.05]);
    ylim=get(gca,'ylim');
    plot([1;1]*fredge',ylim','k:');                         % mark frame boundaries
    hold off
    % v_texthvc(0.02,0.02,sprintf('SNR=%.1f dB',db(sum(s.^2)/sum((sr-s).^2),'p')),'LLk');
    xlabel('Sample Number');
    ylabel('Output Signal + Error');
    title([sprintf('SNR=%.1f dB: ',db(sum(s.^2)/sum((sr-s).^2),'p')) partxt]);
end
v_tilefigs;