% testsynth test with synthetic source signal
%addpath('cadiagit');        % include functions from github
addpath ../utils/
addpath ../tools/

par=projParam;                  % define computer-dependent parameters
pardb=par.db;
maxoff=pardb.maxoff;
% define source of vocal tract parameters
tf = '/TIMIT/TRAIN/DR1/FETB0/SX428.WAV';
tlpc=[0.42 0.56]; % time interval for lpc (a single vowel)
fx=197; % pitch
lpcord=18;
%
% Analyse the TIMIT snippet
%
[stim,fs]=gettimit(tf);
spr=filter([1 -1],1,stim); % differentiate to preemphasise speech
ix=round(tlpc*fs);
ar=v_lpcauto(spr(ix(1):ix(2)),lpcord);
[pf,fpf]=v_lpcar2pf(ar,256);
pfdb=10*log10(pf);
save('lpcar','ar','tf','tlpc','fs');
% figure(11);
% plot(fpf*fs,pfdb);
% xlabel(['Frequency (' v_xticksi 'Hz)']);
% ylabel('LPC filter gain (dB)');
%
% generate some data
%
tt4=0:0.25/fs:2; % 2-second time axis at 4*fx
ud4=v_glotlf(1,tt4*fx); % glottal derivative
ud=resample(ud4,1,4); % subsample to remove aliasing
tt=(0:length(ud)-1)/fs;
s=filter(1,ar,ud(:));
ns=length(s);

%
% derived parameters
%
periodlimk=round(fs./pardb.pitchlim);                     % {max targt mmin} pitch periods (sampgcifracles)
%
% now determine the frame limits
%
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
framemid=[0.5 0.5]*framekk/fs;                      % times of frame centres
%
% now perform transform
%
nfft=2*ceil(periodlimk(1)/2);                       % length of fft (force to be an even number)
[sdft,fax]=gs_stft(s,framex,nfft);
fax=fax*fs;                                         % convert frequency axis to Hz
nfftp=length(fax);                                  % number of positive frequencies
nframe=size(sdft,2);
%
% Use various uniform frames for comparison
%
uwin='Rmn'; % codes for {rectangle, hamming, hanning} windows
uwinl={'Rect','Hamm','Hann'};
% uopt=[1 1; 2 1; 2 2; 2 4]; % Gives [window type, nfft/length]
uopt=[1 1; 3 1; 3 4]; % Gives [window type, nfft/length]
uopt=[3 6; 3 1];
nuopt=size(uopt,1); % number of uniform options
dftus=cell(nuopt,2);
for i=1:nuopt
    dftus{i,1}=v_stftw(s,round(nfft/uopt(i,2)),['t' uwin(uopt(i,1))],1,nfft);
    dftus{i,2}=['STFT ' num2str(round(1000*nfft/(uopt(i,2)*fs))) ' ms ' uwinl{uopt(i,1)}];
end
%
% Now do plots
%
figure(1);
ix=(0:round(3*fs/fx))+round(length(s)/2);
plot(tt(ix),s(ix));
xlabel('Time (s)');
ylabel('Speech signal');
figure(2);
legs=cell(nuopt+2,1);
for i=1:nuopt+2
    switch i
        case 1
            ard=conv(ar,conv([1 -0.95],[1 -0.95])); % apply deemphasis to LPC filter
            vdb=db(v_lpcar2pf(ard,nfftp-1))'/2; 
            legs{1}='Deemph Vocal Tract';   
        case 2
            vdb=db(sdft(1:nfftp,round(nframe/2)));
            legs{2}='Sync STFT';
        otherwise
            ustft=dftus{i-2,1};
            vdb=db(abs(ustft(round(size(ustft,1)/2),:)));        
            legs{i}=dftus{i-2,2};
    end   
    plot(fax,vdb-mean(vdb)-15*i);   
    hold on;
end
hold off
legend(legs);
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%
% v_tilefigs;

