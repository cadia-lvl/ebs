function [melbm,melbu,cfhz,mbm]=stfte2melc(stft,fs,nmel,par)
%  Inputs: stft(nframes,nfft)       complex epoch-based STFT coefficients
%          fs                       sample frequency
%          nmel                     number of mel filters (including DC channel if requested by par.keepDC=1)
%          par                      parameter structure containing optional parameters
%                                       par.fbank=      filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel} ['m']
%                                       par.keepDC      preserve DC as the lowest MEL bin {0, 1} [1]
%                                       par.MELphase    MEL STFT phase calculation: {'true','zero','linear','piecewiselin','piecewiselind','grpdel'} ['piecewiselin']
%                                       par.MELdom      MEL filterbank domain: {'mag', 'pow'} ['pow']
%                                       par.regwt       regularizing factor for phase weights when par.MELphase='piecewiselin' [0.01]
%                                       par.loops       maximum number of iterations when par.MELphase='piecewiselin' [5]
%                                       par.ssqthr      stopping threshold when par.MELphase='piecewiselin' [1.0]
%                                       par.groupdelay  'none'  DFT/DCT control: {'none','dct', ...}
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
% (4) should maybe use meta() input instead of calculating the DFT length directly
% (5) should use the [newish] 'zq' option of v_filtbankm to deal with the par.keepDC option
persistent q0 p2
%
% define default parameters
%
if isempty(q0)
    q0.keepDC=1;                    % preserve DC as the lowest STFT; 0 or 1
    q0.MELphase='piecewiselin';     % MEL STFT phase reconstruction: 'true','zero','linear','piecewiselin'
    q0.MELdom='pow';                % MEL filterbank domain: 'mag', 'pow'
    q0.regwt=0.01;                  % regularizing factor for weights (only for par.MELphase='piecewiselin')
    q0.loops=5;                     % maximum number of iterations (only for par.MELphase='piecewiselin')
    q0.ssqthr=1.00;                 % stopping threshold (only for par.MELphase='piecewiselin')
    q0.fbank='m';                   % filterbank scale: {'b','e','f','m'}={Bark, Erb-rate, Linear, Mel}
    q0.groupdelay='none';            % stft uses DFT rather than DCT
    p2=2*pi;
end
%
% update algorithm parameters
if nargin>=4
    q=v_paramsetch(q0,par);         % update parameters from pp input
else
    q=q0;                           % just copy default parameters
end
if q.keepDC
    fbopt=[q.fbank 'zqdD'];         % filterbank option includes DC output
else
    fbopt=[q.fbank 'dD'];           % filterbank option omits DC output
end

%
% Calculate the mel spectrum
%
nt=size(stft,1);                                                        % nt is number of time frames
if strcmpi(q.groupdelay,'dct') % use DCT instead of DFT
    fbopt=[fbopt 't'];              % add "DCT-input" option for filterbank call
    melbu=[];                       % phase output is null
    nfftprev=-1;                                                            % initialize previous nfft to an impossible value
    for it=1:nt                                                             % loop for all time frames
        nfft=sum(~isnan(stft(it,:)),2);                                     % FFT length for this frame
        newnfft=nfft~=nfftprev;                                             % nfft has changed
        %
        % calculate useful functions of the input stft
        %
        if newnfft                                                          % check if need to recalculate mbm
            %
            % calculate data-independent quantities that depend on nfft
            %
            [mbm,cfhz]=filtbankmd(nmel,nfft,fs,0,fs/2,fbopt);
            mbm=diag(sum(mbm,2).^(-1))*mbm;                             % normalize so that rows sum to 1 to create an interpolation matrix
            cfbin=cfhz*nfft/fs;                                         % MEL centre frequencies in fft bin units (0=dc)
        end
        if strcmpi(q.MELdom,'pow')                                      % check MEL filterbank domain
            melbm(it,:)=sqrt(mbm*abs(stft(it,1:nfft).').^2);                        % MEL filterbank magnitude spectrum in 'pow' dommain
        else
            melbm(it,:)=mbm*stft(it,1:nfft).';                                 % MEL filterbank magnitude spectrum in 'mag' domain
        end
        nfftprev=nfft;                                                  % save nfft to avoid recomputing some quantities if it is unchanged next loop
    end
else
    melbm=zeros(nt,nmel);                                                   % space for mel magnitude spectrum
    if strcmpi(q.MELphase,'true')
        melbu=zeros(nt,size(stft,2));   % space for original phase spectrum
    else
        melbu=zeros(nt,nmel);   % space for mel unwrapped phase spectrum
    end
    nfftprev=-1;                                                            % initialize previous nfft to an impossible value
    for it=1:nt                                                             % loop for all time frames
        nfft=sum(~isnan(stft(it,:)),2);                                     % FFT length for this frame
        newnfft=nfft~=nfftprev;                                             % nfft has changed
        nfftp=1+floor(nfft/2);                                              % number of positive FFT frequency bins
        %
        % calculate useful functions of the input stft
        %
        stftp=stft(it,1:nfftp).';                                           %  input stft for this frame (+ve frequencies only)
        if newnfft                                                          % check if need to recalculate mbm
            %
            % calculate data-independent quantities that depend on nfft
            %
            [mbm,cfhz]=filtbankmd(nmel,nfft,fs,0,fs/2,fbopt);
            mbm=diag(sum(mbm,2).^(-1))*mbm;                             % normalize so that rows sum to 1 to create an interpolation matrix
            cfbin=cfhz*nfft/fs;                                         % MEL centre frequencies in fft bin units (0=dc)
            phlin=lininterps(cfbin,0:nfftp-1,'E');                      % linear interpolation matrix for phases; 'E'=extrapolate ends (nfftp,nmel)
            if q.keepDC                                                 % if mel spectrum includes a DC term ...
                phlin(:,1)=0;                                           % ... don't use DC phase for other bins
                phlin(1,1)=1;                                           % ... except for the DC phase itself
            end
            angs=p2/nfft*(1-nfftp:nfft-nfftp)'*(0:nfftp-1);           % linear phase angles, (nfft,nfftp). The entry angs(k,:) has group delay of nfftp-k samples
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
        switch lower(q.MELphase)
            case 'true'                                                     %%%%% preserve entire true phase spectrum
                melbu(it,:)=angle(stft(it,:));
            case 'zero'                                                     %%%%% set phase spectrum to zero
                if q.keepDC                                                 % if mel spectrum includes a DC term ...
                    melbu(it,1)=angle(stftp(1));                            % make DC phase correct (i.e. 0 or +pi)
                end
            case 'grpdel'                                                   %%%%% interpolate the group delay (note: melbu = group delay in seconds)
                sdfta=angle(stft(it,1:nfft));                               % phase of stft(1,freq)
                sdfta(1)=0;                                                 % zap DC phase
                difa=v_modsym(sdfta-sdfta(previx),p2);                    % phase increment in range +pi
                grpd=zeros(1,nfft);
                grpd(previx)=(difa(previx)+difa)*(nfft/(-4*pi*fs));         % add phase increment for adjacent frequency bins and convert into group delay of seconds
                melbu(it,:)=mbm*grpd(1:nfftp)'; % interpolate group delay to mel bins
            case 'linear'                                                   %%%%% exhaustive search for a linear phase shift
                sdfta=angle(stftp);                                         % phase of stft(1,freq) (+ve frequencies only)
                stftm2=abs(stftp).^2;                                       % stft power for use in the phase-error cost function
                sdftac=cos(sdfta);                                          % cos of true phase angles
                sdftas=sin(sdfta);                                          % sin of true phase angles
                powre=angsc*(stftm2.*sdftac)+angss*(stftm2.*sdftas);        % negative error power (GroupDelay,1)
                [optp,optd]=max(powre);                                     % find optimum group delay in this frame
                melbu(it,:)=2*pi/nfft*cfbin'*(optd-nfftp);                  % calculate linear phases at the MEL centre frequencies to match angs(optd,:) entries
                if q.keepDC                                                 % if mel spectrum includes a DC term ...
                    melbu(it,1)=sdfta(1);                                   % make DC phase correct (i.e. 0 or +pi)
                end
            otherwise                                                       % par.MELphase='piecewiselinu' or'piecewiselin'
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
                % We can also approximate cos(xa-ya)~1-0.5(xa-ya).^2 giving a minimization target (xm.*(xa-ya)).^2          *
                %************************************************************************************************************
                %
                sdfta=angle(stftp);                                         % phase of stft(1,freq) (+ve frequencies only)
                stftm2=abs(stftp).^2;                                       % stft power for use in the phase-error cost function
                mpsqit=sparse(1:nfftp,1:nfftp,sqrt(max(stftm2,max(stftm2)*q.regwt^2)));            % sparse diag matrix square root of regularized phase error weights
                phlinw=mpsqit*phlin;                                    % interpolation matrix weighted by error sensitivity
                ssq=[];
                ssqold=[];
                loopcount=0;
                thresh=q.ssqthr;
                loopmax=q.loops;
                sdftait=unwrap(sdfta); % initialize DFT target unwrapped phase to unwrapped SDFT phase
                if strcmpi(q.MELphase,'piecewiselind') %%%%% piecewise linear unwrapped phase with min-delta regularization
                    phcfit=zeros(nmel,1);                                    % Initialize MEL phases to zero
                    if q.keepDC                                                 % if mel spectrum includes a DC term ...
                        phcfit(1)=sdfta(1);                                   % make DC phase correct (i.e. 0 or +pi)
                    end
                    rega=(q.regwt*max(sum(phlinw,2)))^2;                    % regularization factor
                    phden=phlinw'*phlinw+rega*eye(nmel);                    % denominator matrix ***** note this is only appropriate if nmel<nfftp
                    phnum=phlinw'*mpsqit;                                   % numerator matrix: weighted interpolation
                    % warning('off');                                       % suppress the warnings that happen because phlin is rank-deficient for high nmel
                    phcfit=phden\(phnum*sdftait+rega*phcfit);               % calculate initial phases that minimize weighted error and also minimizes delta
                    while ~numel(ssqold) || ssq<ssqold*thresh && loopcount<loopmax      % loop until the error no longer reduces (should perhaps be the weighted error here)
                        sdftrait=phlin*phcfit;                              % interpolated DFT-resolution unwrapped phase
                        sdftait=v_modsym(sdfta,p2,sdftrait);                % update unwrapped phase of stft(freq,1) to be in range sdftrait+-pi(+ve frequencies only)
                        ssqold=ssq;                                         % remember previous ssq
                        ssq=sum((sdftait-sdftrait).^2);                     % unweighted sum of squared unwrapped phase errors
                        loopcount=loopcount+1;
                    end
                    % warning('on');                                        % re-enable warnings
                    % fprintf(2,'loopcount=%d, ssqreduction=%.1f\n',loopcount,ssq0/ssq);
                else                                                        %%%%% piecewise linear unwrapped phase with smoothness regularization
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % q.MELphase='piecewiselin' finds the least squares solution to the piecewise linear phase reconstruction equation          %
                    % phlin*phcfit=angle(sdft) modulo 2pi. Because of the modulo 2pi term, an iterative scheme is used in which the multiple    %
                    % of 2pi to add to each element of sdfta is updated at each iteration to minimize the squared error. The right-hand-side    %
                    % of the equation is initialized to unwrap(angle(sdft)) for the first iteration.                                            %
                    % The target equation defined above is modified in two respects: (a) each row is multiplied by a weight to reflect its      %
                    % contribution to the phase cost function (see above) and (b) an additional smoothness regularization is incorporated       %
                    % which is needed if any rows of phlin are all-zero (i.e. if any mel bins are unused in the reconstruction).                %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    smcoef=[0.5 -1 0.5]*sqrt(sum(phlinw(:).^2,1)*q.regwt/(nmel-2)); % weighting of smoothness coefficients
                    if q.keepDC                                                     % if preserving DC then ignore DC in regularization
                        amat=[phlinw; toeplitz(zeros(nmel-3,1),[0 smcoef zeros(1,nmel-4)])];
                    else                                                            % ... else include the DC term
                        amat=[phlinw; toeplitz([smcoef(1); zeros(nmel-3,1)],[smcoef zeros(1,nmel-3)])];
                    end
                    % warning('off');                                               % suppress the warnings that happen because phlin is rank-deficient for high nmel
                    phcfit=amat\[mpsqit*sdftait; zeros(nmel-2-q.keepDC,1)];         % calculate initial phases that minimize weighted error to unwrapped SDFT phases
                    while ~numel(ssqold) || ssq<ssqold*thresh && loopcount<loopmax  % loop until the error no longer reduces (should perhaps be the weighted error here)
                        sdftrait=phlin*phcfit;                                      % interpolated full-resolution unwrapped phase
                        sdftait=v_modsym(sdfta,p2,sdftrait);                        % update unwrapped phase of stft(freq,1) to be in range sdftrait+-pi(+ve frequencies only)
                        ssqold=ssq;                                                 % remember previous ssq
                        ssq=sum((sdftait-sdftrait).^2);                             % unweighted sum of squared unwrapped phase errors
                        phcfit=amat\[mpsqit*sdftait; zeros(nmel-2-q.keepDC,1)];     % calculate updated phases that minimize weighted error
                        loopcount=loopcount+1;
                    end
                    sdftraux=phlin*phcfit;                              % interpolated full-resolution unwrapped phase
                    sdftrax=v_modsym(sdftraux,p2);                      % wrapped reconstructed phase in range +-pi
                end
                melbu(it,:)=phcfit;                                     % save updated phases
        end
        nfftprev=nfft;                                                  % save nfft to avoid recomputing some quantities if it is unchanged next loop
    end
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