function [melbm,melbu,cfhz,mbm]=stfte2melc(stft,fs,nmel,par)
%  Inputs: stft(nframes,nfft)       complex epoch-based STFT coefficients
%          fs                       sample frequency
%          nmel                     number of mel filters (including DC channel if requested by par.keepDC=1)
%          par                      parameter structure containing optional parameters
%                                       par.keepDC      preserve DC as the lowest MEL bin {0, 1} [1]
%                                       par.MELphase    MEL STFT phase calculation: {'true','zero','linear','piecewiselin','grpdel'} ['piecewiselin']
%                                       par.MELdom      MEL filterbank domain: {'mag', 'pow'} ['pow']
%                                       par.regwt       regularizing factor for phase weights when par.MELphase='piecewiselin' [0.01]
%                                       par.loops       maximum number of iterations when par.MELphase='piecewiselin' [5]
%                                       par.ssqthr      stopping threshold when par.MELphase='piecewiselin' [1.0]
%
% Outputs: melbm(nframes,nmel)      Positive real-valued Mel filterbank outputs (always in magnitude domain).
%          melbu(nframes,nmel)      Real-valued Mel filterbank unwrapped phase angles
%                                       If par.MELphase='true' then melbu=angle(stft) with dimension (nfft,nframes).
%          cfhz(1,nmel)             Mel-bin centre frequencies in Hz (0=DC)
%          mbm(nmel,nfftp)          Sparse transformation matrix: melbm=mbm*abs(stft(1:nfftp,:)) for 'mag' option
%
%
% Bugs/suggestions:
% (1) should make filterbank peaks always equal 1 or 2 regardless of algorithm
% (2) could conceivably make it cope with non-hermitian stft (from a complex signal)
% (3) par.MELphase='grpdel' might not work with par.keepDC=1 option + doesn't really handle DC right anyway
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.keepDC=1;                    % preserve DC as the lowest STFT; 0 or 1
    q0.MELphase='piecewiselin';     % MEL STFT phase reconstruction: 'true','zero','linear','piecewiselin'
    q0.MELdom='pow';                % MEL filterbank domain: 'mag', 'pow'
    q0.regwt=0.01;                 % regularizing factor for weights (only for par.MELphase='piecewiselin')
    q0.loops=5;                     % maximum number of iterations (only for par.MELphase='piecewiselin')
    q0.ssqthr=1.00;                 % stopping threshold (only for par.MELphase='piecewiselin')
end
%
% update algorithm parameters
if nargin>=4
    q=v_paramsetch(q0,par);             % update parameters from pp input
else
    q=q0;                           % just copy default parameters
end
%
% Calculate the mel spectrum
%
nt=size(stft,1);                                                        % nt is number of time frames
melbm=zeros(nt,nmel);                                                   % space for mel magnitude spectrum
if strcmpi(q.MELphase,'true')
    melbu=zeros(nt,size(stft,2));   % space for original phase spectrum
else
    melbu=zeros(nt,nmel);   % space for mel unwrapped phase spectrum
end
nfftprev=-1;                                                            % initialize previous nfft
for it=1:nt                                                             % loop for all time frames
    nfft=sum(~isnan(stft(it,:)),2);                                     % FFT length for thisframe
    newnfft=nfft~=nfftprev;                                             % nfft has changed
    nfftp=1+floor(nfft/2);                                              % number of positive FFT frequency bins
    %
    % calculate useful functions of the input stft
    %
    stftp=stft(it,1:nfftp).';                                             %  input stft for this frame (+ve frequencies only)   
    if newnfft                                                          % check if need to recalculate mbm
        %
        % calculate data-independent quantities that depend on nfft
        %
        [mbm,cfhz]=v_filtbankm(nmel-q.keepDC,nfft,fs,0,fs/2,'m');
        mbm=diag(sum(mbm,2).^(-1))*mbm;                             % normalize so that rows sum to 1 to create an interpolation matrix
        if q.keepDC                                                 % add extra STFT bin if KEEPdc=1
            mbm(:,1)=0;                                             % eliminate references to DC FFT bin
            mbm=vertcat(sparse(1,1,1,1,nfftp),mbm);                 % preappend an extra row to preserve the DC value
            cfhz=[0 cfhz];                                          % pre-append a centre-frequency for DC
        end
        cfbin=cfhz*nfft/fs;                                             % MEL centre frequencies in fft bin units (0=dc)
        phlin=lininterps(cfbin,0:nfftp-1,'E');              % linear interpolation matrix for phases (nfftp,nmel)
        if q.keepDC                                         % if mel spectrum includes a DC term ...
            phlin(:,1)=0;                                   % ... don't use DC phase for other bins
            phlin(1,1)=1;                                   % ... except for the DC phase itself
        end
        angs=2*pi/nfft*(1-nfftp:nfft-nfftp)'*(0:nfftp-1);           % linear phase angles, (nfft,nfftp). The entry angs(k,:) has group delay of nfftp-k samples
        angsc=cos(angs);                                            % cos of linear phase angles
        angss=sin(angs);                                            % sin of linear phase angles
        previx=[nfft 1:nfft-1];                                     % index of previous element of frame for group delay calculation
    end
    if strcmpi(q.MELdom,'pow')                                      % check MEL filterbank domain
        melbm(it,:)=sqrt(mbm*abs(stftp).^2);                        % MEL filterbank magnitude spectrum in 'pow' dommain
    else
        melbm(it,:)=mbm*abs(stftp);                                 % MEL filterbank magnitude spectrum in 'mag' domain
    end
    %
    % Now calculate phase spectrum
    %
    if strcmpi(q.MELphase,'true')
        melbu(it,:)=angle(stft(it,:));                              % preserve entire true phase spectrum
    elseif strcmpi(q.MELphase,'zero')
        if q.keepDC                                                 % if mel spectrum includes a DC term ...
            melbu(it,1)=angle(stftp(1));                            % make DC phase correct (i.e. 0 or +pi)
        end
    elseif strcmpi(q.MELphase,'grpdel')                             % output group delay in seconds
        sdfta=angle(stft(it,1:nfft));                               % phase of stft(1,freq)
        sdfta(1)=0; % zap DC phase
        difa=v_modsym(sdfta-sdfta(previx),2*pi);                    % phase increment in range +-pi
        grpd=zeros(1,nfft);
        grpd(previx)=(difa(previx)+difa)*(nfft/(-4*pi*fs));         % add phase increment for adjacent frequency bins and convert into group delay of seconds
        melbu(it,:)=mbm*grpd(1:nfftp)';
    else
        %************************************************************************************************************
        % Phase cost function:                                                                                      *
        %                                                                                                           *
        % For one frame, if true stft is x = xm.*exp(1i*xa) and estimated stft is y = ym.*exp(1i*ya)                *
        % the the estimation error is                                                                               *
        %    sum(|x-y|.^2) = (x-y)'*(x-y)                                                                           *
        %            = (xm exp(1i xa) - ym exp(1i ya))'*(xm exp(1i xa) - ym exp(1i ya))                             *
        %            = xm'*xm + ym'*ym - 2*(xm.*ym)'*cos(xa-ya)                                                     *
        %                                                                                                           *
        % So, if the output magnitude is fixed, we minimize the error by choosing ya to maximize                    *
        %     (xm.*ym)'*cos(xa-ya) = (xm.*ym)'*(cos(xa).*cos(ya)+sin(xa).*sin(ya))                                  *
        %                          = (xm.*ym.*cos(xa))'*cos(ya)+(xm.*ym.*sin(xa))'*sin(ya)                          *
        % We approximate xm.*ym by xm.^2 to avoid the phase dependency on the magnitude reconstruction algorithm    *
        %************************************************************************************************************
        %
        % first do an exhaustive search for a linear phase solution
        %
        sdfta=angle(stftp);                                                 % phase of stft(1,freq) (+ve frequencies only)
        stftm2=abs(stftp).^2;                                   % stft power for use in the phase-error cost function
        sdftac=cos(sdfta);                                      % cos of true phase angles
        sdftas=sin(sdfta);                                      % sin of true phase angles
        powre=angsc*(stftm2.*sdftac)+angss*(stftm2.*sdftas);    % negative error power (GroupDelay,1)
        [optp,optd]=max(powre);                                 % find optimum group delay in this frame
        % angs=2*pi/nfft*(1-nfftp:nfft-nfftp)'*(0:nfftp-1);
        melbu(it,:)=2*pi/nfft*cfbin'*(optd-nfftp);              % calculate linear phases at the MEL centre frequencies to match angs(optd,:) entries
        if q.keepDC                                             % if mel spectrum includes a DC term ...
            melbu(it,1)=sdfta(1);                               % make DC phase correct (i.e. 0 or +pi)
        end
        loopmax=q.loops;
        if strcmpi(q.MELphase,'piecewiselin')
            % mpsqit=sparse(1:nfftp,1:nfftp,sqrt((stftm2+q.regwt*max(stftm2)))); % diag matrix regularized square root of phase error weights
            mpsqit=sparse(1:nfftp,1:nfftp,sqrt(stftm2)); % diag matrix regularized square root of phase error weights
            phcfit=melbu(it,:)';                             % Initialize MEL weights to linear phase values
            phlinw=mpsqit*phlin; % interpolation matrix weighted by error sensitivity
            rega=(q.regwt*max(sum(phlinw,2)))^2; % regularization factor
            phden=phlinw'*phlinw+rega*eye(nmel);
            phnum=phlinw'*mpsqit;
            ssq=[];
            ssqold=[];
            % warning('off');                                 % suppress the warnings that happen because phlin is rank-deficient for high nmel
            loopcount=0;
            thresh=q.ssqthr;
            while ~numel(ssqold) || ssq<ssqold*thresh && loopcount<loopmax            % loop until the error no longer reduces (should perhaps be the weighted error here)
                sdftrait=phlin*phcfit;                      % interpolated full-resolution unwrapped phase
                sdftait=mod(sdfta-sdftrait+pi,2*pi)+sdftrait-pi;  % unwrapped phase of stft(freq,1) (+ve frequencies only)
                ssqold=ssq;                                 % remember previous ssq
                ssq=sum((sdftait-sdftrait).^2);             % unweighted sum of squared unwrapped phase errors
                phcfit=phden\(phnum*sdftait+rega*phcfit);         % calculate updated phases that minimize weighted error
                % if loopcount==0
                %     ssq0=ssq;
                % end
                loopcount=loopcount+1;
            end
            % warning('on');                                  % re-enable warnings
            melbu(it,:)=phcfit;                             % save updated phases
            % fprintf(2,'loopcount=%d, ssqreduction=%.1f\n',loopcount,ssq0/ssq);
        end
    end
    nfftprev=nfft;                                              % save nfft to avoid recomputing some quantities if it is unchanged next loop
end
if ~nargout  
    subplot(2,1,2);
    imagesc(melbu'*(180/pi));
    axis 'xy';
    colorbar;
    xlabel('Frame Num');
    ylabel('Mel bin');
    cblabel('U-phase (deg)')
    subplot(2,1,1);
    melbmmax=max(melbm(:));
    melbmdb=20*log10(max(melbm,1e-3*melbmmax));
    imagesc(melbmdb');
    axis 'xy';
    colorbar;
    ylabel('Mel bin');
    cblabel('Mag (dB)');
end