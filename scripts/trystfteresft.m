% trial of stfteresft
close all;
par=stftepinit(); % initialize all parameters to default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation parameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.interpdom=  'magcph';       % Interpolatione domain: {'cplx','magcph','crmcph'}
par.interpseq=  1;              % 0=interpolate T and F jointly, 1=interpolate T and F sequentially with an intermediate complex spectrum
par.interpp=    5;              % length of intermediate filter
%
timit='D:/OneDrive - Imperial College London/work/data/speech/timit/timit/';    % path to timit sub-folder of timit CD
tf='TRAIN/DR6/MSJK0/SX246.WAV';         % filename=TIMIT file to use
[s,fs,wrd,phn]=gettimit(tf,'n',timit);  % read the file
ns=length(s);                           % number of speech samples
[fein,gci]=gs_frames(s,fs,par);         % index of last sample in each frame
gcik=round(gci*fs+1);                   % gci positions (sample numbers)
frameend=smoothframes(s,fein,par);      % smooth the frame pitch contour
nframe=length(frameend);
framestart=[1 frameend(1:end-1)+1];     % frame starting sample numbers
framecent=0.5*(framestart+frameend);    % frame centres (fractional sample numbers)
metain=[framestart(:) frameend(:)-framestart(:)+1 repmat([0 1 0 0 0 0],nframe,1)]; % initialize metadata
[stftx,metax,gdsh,grpdx,dgrpdx]=stfte(s,metain); % do the stft
% now do the interpolation
%
nfft=round(fs/100);             % length of interpolated DFT [100 Hz]
flen=round(fs*6e-3);            % length of output frames in samples [6 ms]
[stftv,metav]=stfteresft(stftx,metax,par,nfft,flen); % interpolate onto a fixed grid
% for comparison use stftegridz
[stftg,metag]=stftegridz(stftx,metax,metav(:,1:3),par);
%
% now do plots
%
axlinkxt=[];                                            % time axis linking list
axlinkxf=[];                                            % frequency axis linking list
axlinkxytf=[];                                          % time-frequency axis linking list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: time-domain waveform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);                                              
plot((1:length(s))/fs,s);
v_axisenlarge([-1 -1.05]);
axlinkxt=[axlinkxt gca]; % link time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2: frame lengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
plotstftemeta('lT',metax,fs); % plot frame lengths
axlinkxt=[axlinkxt gca]; % link time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3: spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
v_spgrambw(s,fs,'pJcwaAtTi',[],[],[],[],{wrd phn}); % plot spectrogram
axlinkxt=[axlinkxt gca]; % link time axis
axlinkxytf=[axlinkxytf gca]; % link time-frequency axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4: Original Epoch-based Mag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
stfteplot(stftx,metax,fs,40);
title('Original Epoch-based Mag');
axlinkxytf=[axlinkxytf gca]; % link time-frequency axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 8:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8);
stfteplot(stftv,metav,fs,40);
title(sprintf('New TF-interpolated (%.0f ms frames)',1000*flen/fs));
axlinkxytf=[axlinkxytf gca]; % link time-frequency axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 7:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
stfteplot(stftg,metag,fs,40);
title(sprintf('Old TF-interpolated (%.0f ms frames)',1000*flen/fs));
axlinkxytf=[axlinkxytf gca]; % link time-frequency axis
%
v_tilefigs;
if ~isempty(axlinkxt)
    linkaxes(axlinkxt,'x');
end
if ~isempty(axlinkxf)
    linkaxes(axlinkxf,'x');
end
if ~isempty(axlinkxytf)
    linkaxes(axlinkxytf,'xy');
end