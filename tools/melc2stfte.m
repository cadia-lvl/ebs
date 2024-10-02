function stft=melc2stfte(melbm,melbu,fs,meta,par)
%
%  Inputs: melbm(nframes,nmel)     Positive real-valued Mel filterbank outputs (always in magnitude domain).
%          melbu(nframes,nmel)     Real-valued Mel filterbank unwrapped phase angles
%                                       If par.MELphase='true' then melbu=angle(stft) with dimension (nfft,nframes).
%          cfbin(nmel)               Mel-bin centre frequencies in units of fft bins (0=DC)
%          fs                       Sample frequency
%           meta(nframes,6)         We only use meta(:,3) which is the fft length
%          par                      parameter structure containing optional parameters
%                                       par.fbank=      filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel} ['m']
%                                       par.keepDC      preserve DC as the lowest MEL bin {0, 1} [1]
%                                       par.invmethod   MEL inversion method: {'pinv','transp','linterp','fbxi'} ['linterp']
%                                       par.MELphase    MEL STFT phase calculation: {'true','zero','linear','piecewiselin'} ['piecewiselin']
%                                       par.MELdom      MEL filterbank domain: {'mag', 'pow'} ['pow']
%
% Outputs: stft(nframes,nfftmax)       complex epoch-based STFT coefficients
%
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.keepDC=1;                    % preserve DC as the lowest STFT; 0 or 1
    q0.invmethod='linterp';         % Magnitude reconstruction method: 'pinv' or 'linterp'
    q0.MELphase='piecewiselin';     % MEL STFT phase reconstruction: 'true','zero','linear','piecewiselin'
    q0.MELdom='pow';                % MEL filterbank domain: 'mag', 'pow'
    q0.fbank='m';                   % filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel}
end
%
% update algorithm parameters
if nargin>=5
    q=v_paramsetch(q0,par);             % update parameters from pp input
else
    q=q0;                           % just copy default parameters
end
if q.keepDC
    fbopt=[q.fbank 'zqdD']; % filterbank option includes DC output
else
    fbopt=[q.fbank 'dD'];  % filterbank option omits DC output
end
%
% sort out input parameters
%
[nt,nmel]=size(melbm);
fftlen=meta(:,3);
nfftmax=max(fftlen);               % ... allow room for biggest fft
%
% Calculate inverse mel spectrum
%
stft=NaN(nt,nfftmax);           % initialize stft to all NaN
nfftprev=-1;                                                            % initialize previous nfft
for it=1:nt                                                             % loop for all time frames
    nfft=fftlen(it);                                              % fft length for this frame
    newnfft=nfft~=nfftprev;                                             % nfft has changed
    nfftp=1+floor(nfft/2);                                              % number of positive FFT frequency bins
    melbmi=melbm(it,:)';
    melbui=melbu(it,:)';
    if newnfft
        %
        % here goes data-independent code that depends on fft length
        %

        %
        % create the inverse mbm matrix
        %
        switch q.invmethod                                              % pseudo-inverse version
            case 'pinv'
                [mbm,cfhz]=v_filtbankm(nmel,nfft,fs,0,fs/2,fbopt);
                mbm=diag(sum(mbm,2).^(-1))*mbm;                         % normalize so that rows sum to 1 to create an interpolation matrix
                % if q.keepDC                                             % add extra STFT bin if KEEPdc=1
                %     mbm(:,1)=0;                                         % eliminate references to DC FFT bin (set column 1 to zero)
                %     mbm=vertcat(sparse(1,1,1,1,nfftp),mbm);             % preappend an extra row to preserve the DC value
                % end
                imbm=pinv(full(mbm));
            case 'fbxi'                                                 % use inverse from v_filtbankm
                [mbm,cfhz,imbm]=v_filtbankm(nmel,nfft,fs,0,fs/2,fbopt);
            case 'transp'                                               % transpose version
                [mbm,cfhz]=v_filtbankm(nmel,nfft,fs,0,fs/2,fbopt);
                imbm=diag(sum(mbm,1).^(-1))*mbm';                       % normalize so that rows sum to 1 to create an interpolation matrix
            case 'linterp'                                              %  linear interpolation version
                nmelx=nmel-q.keepDC;
                switch q.fbank
                    case 'm'
                        cfhz=v_mel2frq((0:nmelx+1)*v_frq2mel(fs/2)/(nmelx+1));
                    case 'b'
                        cfhz=v_bark2frq((0:nmelx+1)*v_frq2bark(fs/2)/(nmelx+1));
                    case 'e'
                        cfhz=v_erb2frq((0:nmelx+1)*v_frq2erb(fs/2)/(nmelx+1));
                    case 'f'
                        cfhz=(0:nmelx+1)*0.5*fs/(nmelx+1);
                end
                cfhz=cfhz(2-q.keepDC:nmelx+1); % remove DC bin if q.keepDC==0
                cfbin=cfhz*nfft/fs; % centre frequencies in fft bin units (0=dc)
                imbm=lininterps(cfbin,0:nfftp-1);
                if q.keepDC                                             % if mel spectrum includes a DC term ...
                    imbm(1,:)=0;                                        % ... eliminate other influences on DC bin
                    imbm(:,1)=0;                                        % ... don't use DC mel-bin for any other fft bin
                    imbm(1,1)=1;                                        % ... map DC bin directly without scaling
                end
            otherwise
                error(['filterbank inversion method ' q.invmethod ' not defined'])
        end
        cfbin=cfhz*nfft/fs;                                             % MEL centre frequencies in fft bin units (0=dc)
    end
    if strcmpi(q.MELdom,'pow')                                   % check MEL filterbank domain
        stftp=sqrt(max(imbm*melbmi.^2,0));                           % reconstructed magnitude stft in 'pow' domain
    else
        stftp=max(imbm*melbmi,0);                                    % reconstructed magnitude stft in 'mag' domain
    end
    %
    % now interpolate the phase
    %
    if strcmpi(q.MELphase,'true')
        stftp=stftp.*exp(1i*melbui(1:nfftp));                                 % reconstruct with true phases
    elseif strcmpi(q.MELphase,'zero')
        if q.keepDC                                                 % if mel spectrum includes a DC term ...
            stftp(1)=stftp(1).*(2*(melbui(1)==0)-1);           % correct the sign of the DC term
        end
    else                                                            % linearly interpolate phases except of the DC term
        phlin=lininterps(cfbin,0:nfftp-1,'E');                      % linear interpolation matrix for phases (nfftp,nmel)
        if q.keepDC                                                 % if mel spectrum includes a DC term ...
            phlin(:,1)=0;                                           % ... don't use DC phase for other bins
            phlin(1,1)=1;                                           % ... except for the DC phase itself
        end
        stftp=stftp.*exp(phlin*(1i*melbui));                         % interpolated full-resolution phase
    end
    stft(it,1:nfft)=[stftp; conj(stftp(1+nfft-nfftp:-1:2,:))];                 % reconstruct negative frequencies
    nfftprev=nfft;                                              % save nfft to avoid recomputing some quantities if it is unchanged next loop
end
if ~nargout
    rgb1=v_colormap('v_thermliny','k');                     % get colormap for magnitude
    rgb2=v_colormap('v_circrby','k');                       % get colormap for phase
    nfftp=1+floor(nfftmax/2);                                % number of positive frequencies from largest fft
    faxz=(0:nfftp-1)/nfftmax;                                % frequency axis for spectrogram
    astft=zeros(nt,nfftp);                      % space for the interpolated stft
    for i=1:nt
        astft(i,:)=stft(i,1+round((0:nfftp-1)*meta(i,3)/nfftmax));      % interpolate to fill astft by replicating
    end
    minval=max(abs(astft(:)))*1e-2;
    subplot(2,1,2);
    imagesc(1:nt,faxz,mod(angle(astft)'*180/pi+360/128,360)-360/128,360*[-1 127]/128);    % ensure data range is 360*[-1 127]/128
    colorbar;
    cblabel('Phase (Deg)');
    colormap(gca,rgb2);
    axis('xy');
    xlabel('Frame Num');
    ylabel('Freq\div{}fs');
    subplot(2,1,1);
    imagesc(1:nt,faxz,db(max(abs(astft)',minval)));
    colorbar;
    axis('xy');
    ylabel('Freq\div{}fs');
    cblabel('Mag (dB)');
    colormap(gca,rgb1);
end
