function [stft,meta,gdsh,grpd,dgrpd]=stfte(s,metain,maxfft,par)
% Epoch-based stft
%
%  Inputs: s(ns,1)                  real-valued signal
%          metain(nframe,2)         input metadata: meta(*,:)=[first-sample,frame-length]
%          maxfft                   max size of frame in stft domain [default: max(metain(:,2))]
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.window      'r'         Window: {'r','n','m'} = {Rectangular, Hanning, Hamming}
%                                       par.offset      'none'      offset removal: {'none','mean','ends'}
%                                       par.scale       'none'      scaling method: {'none','peakabs','rms'}
%                                       par.pad         'none'      zero-padding method: {'none','zero','ends'}
%                                       par.gcifrac     0           position of GCI in analysis frame. Only used if the par.groupdelay='gcif' option is selected.
%                                       par.groupdelay  'none'      linear phase component: {'none','dct','ewgd','cplx','phgr','gcif','gpdf','fmnb','xcor'} + optional 'int' suffix
%                                       par.fmbound     [0.3 0.5]   group delay bounds for par.groupdelay='fmnb'
%
% Outputs: stft(nframe,maxfft)      complex STFT coefficients
%          meta(nframe,10)           output metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples), EWGD (samples), EWPD (samples), DC phase (rad), DC group delay (samples)]
%                                   Note that if there are less than 3 output arguments, meta(:,10) will be set to zero.
%          grpd(nframe,maxfft)      group delay in samples. This is calculated from the slope of a quadratic fitted to three consecutive points along the frequency axis.
%          dgrpd(nframe,maxfft)     derivative of group delay in samples^2. This is calculated from a quadratic fitted to three consecutive points along the frequency axis.
%
%
% Note that the stft phase (modulo 2pi) for frame i can be recovered exactly from meta(i,:) and either gdsh(i,:) or dgrpd(i,:) as:
%           nfft=meta(i,3);
%           ddc=-2*pi*meta(i,10)/nfft;
%           phase1=cumsum([meta(i,9) gdsh(i,1:nfft-1)*-2*pi./nfft],2);
%           phase2=cumsum(cumsum([ddc -(2*pi./nfft)^2*dgrpd(i,1:nfft-1)],2),2)-ddc+meta(i,9);
%
% Transformations are applied in the order window, offset, scale, pad, dft, groupdelay. If pad option is 'ends', the group delay can exceed the length of the unpadded frame.
%
% Algorithm Options:
%
%   par.window      'r'         Rectangular [default]
%                   'n'         Hanning window applied before offset/scale/pad operations
%                   'm'         Hamming window applied before offset/scale/pad operations
%   par.windowmode  'E'         window mode (see v_windows.m) [currently ignored for the rectangular window]
%   par.offset      'none'      No offset-removal is performed so meta(:,4)=0 [default]
%                   'mean'      The mean of each frame is subtracted from the frame
%                   'ends'      The average of the first and last samples in each frame is subtracted from the frame
%   par.scale       'none'      No scaling is performed so meta(:,5)=1 [default]
%                   'peakabs'   Scale each frame so the peak absolute value equals 1
%                   'rms'       Scale each frame so that the root-mean-square value equals 1
%                   'len'       Scale each frame by its unpadded length (total spectral energy = average energy per sample)
%                   'sqlen'     Scale each frame by the square root of its unpadded length (average spectral energy = average energy per sample)
%   par.pad         'none'      No padding is performed so dft-length for each frame equals the frame length and meta(:,3)=meta(:,2) [default]
%                   'zero'      Pad each frame with zeros to a length of maxfft
%                   'ends'      Padd each frame with the average of the end values to a length of maxfft
%   par.gcifrac                 position of GCI in analysis frame used only by par.groupdela='gcif' option below [0]
%   par.groupdelay  'none'      No group delay compensation is performed so meta(:,6)=0 (i.e. zero samples delay) [default]
%                   'dct'       Use DCT rather than DFT (so no phases)
%                   'ewgd'      Take group delay equal to the centre of gravity of the energy of the frame after padding
%                   'cplx'      As ewgd but convert position in frame to a complex number before doing the energy-weighted
%                               average. This mirrors the circularity of the DFT.
%                   'phgr'      Calculate the enegrgy-weighted phase decrement beween successive frequency bins in the DFT output
%                   'gcif'      Take group delay equal to par.gcifrac multiplied by the frame length, meta(:,2). Note that par.gcifrac defaults to 0 unless it is explicitly specified.
%                   'gpdf'      Take group delay equal to par.gpdfrac multiplied by the frame length, meta(:,2) [default=0.3].
%                   'fmnb'      Find optimum group delay subject to bounds par.fmbound as fraction of frame length [default bounds = [0.3 0.5]].
%                               'fmnb' minimizes an energy-weighted average of 1-cos(phi). It is quite slow.
%                   'xcor'      for fixed frame size only; maximize the cross-correlation between successive frames
%                   '****int'   As above but rounded to an integer number of samples where '****' is one of the previous options
%   par.fmbound                 Group delay bounds as fraction of frame length when using par.groupdelay='fmnb' option [default = [0.3 0.5]]
%
% Bugs/Suggestions:
% (1) Group delay calculation could possibly include the option of a phase shift of pi in consecutive
%     frequency bins due to magnitude sign change (especially at DC and Nyquist frequency).
%
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.offset='none';                   %  offset removal: {'none','mean','ends'}
    q0.scale='none';                    %  scaling method: {'none','peakabs','rms','len','sqlen'}
    q0.pad='none';                      % zero-padding method: {'none','zero','ends'}
    q0.groupdelay='none';               % linear phase component: {'none','dct','ewgd','cplx','phgr','gcif','gpdf','fmnb' + optional 'int' suffix}
    q0.gcifrac=0;                       % fraction of frame length for par.groupdelay='gpdf' option
    q0.gpdfrac=0.3;                     % fraction of frame length for par.groupdelay='gpdf' option
    q0.fmbound=[0.3 0.5];               % bounds for group delay option par.groupdelay='fmbd' as fraction of frame length
    q0.window='r';                      % window: 'r'=rectangular, 'n'=hanning, 'm'=hamming
    q0.windowmode='Es';                  % window mode (see v_windows.m)
end
%
% update algorithm parameters
if nargin<4
    q=q0;                               % just copy default parameters
else
    q=v_paramsetch(q0,par);             % update parameters from par input
end
%
nframe=size(metain,1);                  % number of frames
meta=zeros(nframe,10);                   % space for metadata
meta(:,1:2)=metain(:,1:2);              % copy frame start and frame length information
framelens=meta(:,2);                    % length of each frame in samples
if nargin<3 || isempty(maxfft)
    maxfft=max(framelens);              % set maxfft to maximum frame length
else
    framelens=min(framelens,maxfft);    % no framelength can exceed maxfft
end
sfr=zeros(nframe,maxfft);                                       % space for frames (one per row)
stft=NaN(nframe,maxfft);                                        % space for stft; unused entries will be left as NaN
grpd=NaN(nframe,maxfft);                                        % space for group delay; unused entries will be left as NaN
gsh=NaN(nframe,maxfft);                                        % space for shifted group delay; unused entries will be left as NaN
dgrpd=NaN(nframe,maxfft);                                       % space for derivative of group delay; unused entries will be left as NaN
sfrix=repmat(0:maxfft-1,nframe,1)+repmat(meta(:,1),1,maxfft);   % index into s for frame samples
sfrmk=repmat(1:maxfft,nframe,1)<=repmat(meta(:,2),1,maxfft);    % mask for valid values
sfr(sfrmk)=s(sfrix(sfrmk));                                     % enframe the data
if q.window~='r'
    for i=1:nframe % loop through frames (could vectorize for simple windows)
        sfr(i,1:meta(i,2))=sfr(i,1:meta(i,2)).*v_windows(q.window,meta(i,2),q.windowmode);
    end
end
switch q.offset
    case 'mean'
        meta(:,4)=sum(sfr,2)./meta(:,2);                                % find mean of each frame
        sfrrep=repmat(meta(:,4),1,maxfft);
        sfr(sfrmk)=sfr(sfrmk)-sfrrep(sfrmk);                            % subtract mean from valid samples
    case 'ends'
        meta(:,4)=0.5*(sfr(:,1)+sfr((meta(:,2)-1)*nframe+(1:nframe)')); % find average of the end-samples of each frame
        sfrrep=repmat(meta(:,4),1,maxfft);
        sfr(sfrmk)=sfr(sfrmk)-sfrrep(sfrmk);                            % subtract mean from valid samples
end
if strcmp(q.scale,'rms')
    meta(:,5)=sqrt(sum(sfr.^2,2)./meta(:,2));               % find rms in each frame
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
elseif strcmp(q.scale,'peakabs')
    meta(:,5)=max(abs(sfr),[],2);                           % find max of abs(sfr)) in each frame
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
elseif strcmp(q.scale,'sqlen')
    meta(:,5)=sqrt(meta(:,2));                              % scale by the square root of the frame length
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
elseif strcmp(q.scale,'len')
    meta(:,5)=meta(:,2);                                    % scale by the frame length
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
else
    meta(:,5)=1;                                            % default scale factor is unity
end
if strcmp(q.pad,'none')                                     % if no padding
    meta(:,3)=meta(:,2);                                    % dft length equals the frame length
else                                                        % need to include padding
    meta(:,3)=maxfft;                                       % dft length equals maxfft
end
if all(framelens==framelens(1))                             % all frames are the same length so we can do them all at once
    framelen=framelens(1);                                  % constant frame length for all frames
    if strcmp(q.pad,'ends') && framelen~=maxfft             % we need to pad frames with the average of the endpoints
        sfr(:,framelen+1:maxfft)=repmat(sfr(:,[1 framelen])*[0.5;0.5],1,maxfft-framelen);   % pad with average of endpoints
    end
    nfft=meta(1,3);                                         % DFT length of all frames
    if strcmp(q.groupdelay,'dct')                          % do DCT rather than DFT
        stft(:,1:nfft)=v_rdct(sfr(:,1:nfft)')';             % Perform DCT on all frames of sfr (rows)
    else
        stft(:,1:nfft)=fft(sfr(:,1:nfft),nfft,2);           % Perform DFT on all frames of sfr (rows)
        if ~strcmp(q.groupdelay,'none')
            switch q.groupdelay(1:4)
                case 'xcor'                                 % empty for now
                    [ddz,jxz]=max(real(ifft(conj(stft(1:nframe-1,:)).*stft(2:nframe,:),nfft,2)),[],2); % find peak correlation. jxz=1 for zero lag, jxz(i)=2 if frame i+1 lags frame i by 1 sample
                    meta(2:end,6)=mod(cumsum(jxz-1),nfft); % group delay is the cumulative sum of the inter-frame lags
                case 'ewgd'
                    meta(:,6)=sfr(:,1:nfft).^2*(0:nfft-1)'./sum(sfr(:,1:nfft).^2,2); % calculate EWGD for all frames
                case 'cplx'
                    meta(:,6)=mod(angle(sfr(:,1:nfft).^2*exp(2i*pi*(0:nfft-1)'/nfft))*nfft/(2*pi),nfft); % calculate complex group delay for this frame
                case 'phgr'
                    meta(:,6)=angle(sum(stft(:,1:nfft-1).*conj(stft(:,2:nfft)),2)./sum(abs(stft(:,1:nfft-1).*stft(:,2:nfft)),2))*nfft/(2*pi);
                case 'gcif'
                    meta(:,6)=meta(:,2)*q.gcifrac;
                case 'gpdf'
                    meta(:,6)=meta(:,2)*q.gpdfrac;
                case 'fmnb'
                    for i=1:nframe
                        meta(i,6)=fminbnd(@(g) gderr(g,stft(i,1:nfft)),q.fmbound(1)*(framelen-1),q.fmbound(2)*(framelen-1)); % find optimum
                    end
                otherwise
                    error(sprintf('par.goupdelay equals unknown value: %s',q.groupdelay));
            end
            if length(q.groupdelay)>4
                meta(:,6)=round(meta(:,6));                     % round EWGD to an integer
            end
            stft(:,1:nfft)=stft(:,1:nfft).*exp(2i*pi/nfft*meta(:,6)*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1]); % apply group delay (except to Nyquist frequency)
        end
        meta(:,7)=v_modsym(sum(sfr(:,1:nfft).^2.*repmat(0:nfft-1,nframe,1),2)./sum(sfr(:,1:nfft).^2,2),nfft); % calculate energy-weighted group delay from time-domain signal
        meta(:,8)=v_modsym(angle(sum(sfr(:,1:nfft).^2.*repmat(exp(2i*pi*(0:nfft-1)/nfft),nframe,1),2))*nfft/(2*pi),nfft); % calculate energy-weighted phase delay from time-domain signal
        meta(:,9)=pi*(real(stft(:,1))<0); % phase at DC (always 0 or pi)
        if nargout>2                                            % need to calculate group delay also
            stfta=angle(stft(:,1:nfft));                        % calculate phases
            previx=[nfft 1:nfft-1];                             % index of previous element of frame
            difa=v_modsym(stfta-stfta(:,previx),2*pi);          % phase increment in range +-pi
            gdsh(:,previx)=difa*(framelen/(-2*pi)); % group delay in samples at shifted frequencies
            grpd(:,previx)=(difa(:,previx)+difa)*(framelen/(-4*pi));    % add phase increment for adjacent frequency bins and convert into group delay of samples
            dgrpd(:,previx)=v_modsym(difa-difa(:,previx),2*pi)*((-0.25/pi^2)*framelen^2); % derivative of group delay
            meta(:,10)=difa(:,1)*((-0.5/pi)*framelen);
        end
    end
else                                                        % we must process frames individually 'cos varying lengths
    for i=1:nframe
        framelen=framelens(i);                              % length of this frame in samples
        if strcmp(q.pad,'ends') && framelen~=maxfft
            sfr(i,framelen+1:maxfft)=repmat(sfr(i,[1 framelen])*[0.5;0.5],1,maxfft-framelen); % pad with average of endpoints
        end
        nfft=meta(i,3);                                     % DFT length
        if strcmp(q.groupdelay,'dct')                          % do DCT rather than DFT
            stft(i,1:nfft)=v_rdct(sfr(i,1:nfft)')';             % Perform DCT on current frame of sfr (row)
        else
            stft(i,1:nfft)=fft(sfr(i,1:nfft));                  % Perform DFT
            if ~strcmp(q.groupdelay,'none')
                switch q.groupdelay(1:4)
                    case 'xcor'                                 % empty for now
                        error('xcor group delay option not yet implemented for variable-length frames');
                    case 'ewgd'
                        meta(i,6)=sfr(i,1:nfft).^2*(0:nfft-1)'/sum(sfr(i,1:nfft).^2); % calculate EWGD for this frame
                    case 'cplx'
                        meta(i,6)=mod(angle(sfr(i,1:nfft).^2*exp(2i*pi*(0:nfft-1)'/nfft))*nfft/(2*pi),nfft); % calculate complex group delay for this frame
                    case 'phgr'
                        meta(i,6)=mod(angle(stft(i,2:nfft-1)*stft(i,3:nfft)'/sum(abs(stft(i,2:nfft-1).*stft(i,3:nfft))))*nfft/(2*pi),nfft); % omit DC from calculation
                    case 'gcif'
                        meta(i,6)=meta(i,2)*q.gcifrac;
                    case 'gpdf'
                        meta(i,6)=meta(i,2)*q.gpdfrac;
                    case 'fmnb'
                        % meta(i,6)=fminbnd(@(g) gderr(g,stft(i,1:nfft)),q.fmbound(1)*(framelen-1),q.fmbound(2)*(framelen-1)); % find optimum
                        [xx,yy,bb]=fftgdopterr(stft(i,1:nfft));                             % initial call computes fixed quantities for optimization
                        meta(i,6)= fminbnd(@(g) fftgdopterr(g,xx,yy,bb),q.fmbound(1)*(framelen-1),q.fmbound(2)*(framelen-1)); % find optimum
                    otherwise
                        error(sprintf('par.goupdelay equals unknown value: %s',q.groupdelay));
                end                                             % if ~strcmp(q.groupdelay,'none') ... end
                if length(q.groupdelay)>4
                    meta(i,6)=round(meta(i,6));                 % round group delay to an integer
                end
                stft(i,1:nfft)=stft(i,1:nfft).*exp(2i*pi/nfft*meta(i,6)*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1]); % apply group delay (except to Nyquist frequency)
            end
            stfta=angle(stft(i,1:nfft));                        % calculate phases
            previx=[nfft 1:nfft-1];                             % index of previous element of frame
            difa=v_modsym(stfta-stfta(previx),2*pi);            % phase increment in range +-pi
            gdsh(i,previx)=difa*(framelen/(-2*pi)); % group delay in samples at shifted frequencies
            grpd(i,previx)=(difa(previx)+difa)*(framelen/(-4*pi));      % add phase increment for adjacent frequency bins and convert into group delay of samples
            dgrpd(i,previx)=v_modsym(difa-difa(previx),2*pi)*((-0.25/pi^2)*framelen^2); % derivative of group delay
            meta(i,7)=v_modsym(sum(sfr(i,1:nfft).^2.*(0:nfft-1),2)./sum(sfr(i,1:nfft).^2,2),nfft); % calculate energy-weighted group delay from time-domain signal
            meta(i,8)=v_modsym(angle(sum(sfr(i,1:nfft).^2.*exp(2i*pi*(0:nfft-1)/nfft),2))*nfft/(2*pi),nfft); % calculate energy-weighted phase delay from time-domain signal
            meta(i,9)=stfta(1);                                 % phase at DC (always 0 or pi)
            meta(i,10)=difa(1)*((-0.5/pi)*framelen);             % group delay at DC (might not equal grpd(i,1) if stfta(1)=+-pi)
        end                                                     % if strcmp(q.groupdelay,'dct')   ... else ... end
    end                                                         % for i=1:nframe ... end
end                                                             % if all(framelens==framelens(1))  ... else ... end
if ~nargout
    rgb1=v_colormap('v_thermliny','k');                     % get colormap for magnitude
    rgb2=v_colormap('v_circrby','k');                       % get colormap for phase
    nfftp=1+floor(maxfft/2);                                % number of positive frequencies
    faxz=(0:nfftp-1)/maxfft;                                % frequency axis for spectrogram
    astft=zeros(nframe,nfftp);                              % space for the interpolated stft
    for i=1:nframe
        astft(i,:)=stft(i,1+round((0:nfftp-1)*meta(i,3)/maxfft));      % interpolate to fill astft by replicating
    end
    minval=max(abs(astft(:)))*1e-2;
    subplot(2,1,2);
    imagesc(1:nframe,faxz,mod(angle(astft)'*180/pi+360/128,360)-360/128,360*[-1 127]/128);    % ensure data range is 360*[-1 127]/128
    colorbar;
    v_cblabel('Phase (Deg)');
    colormap(gca,rgb2);
    axis('xy');
    xlabel('Frame Num');
    ylabel('Freq\div{}fs');
    subplot(2,1,1);
    imagesc(1:nframe,faxz,db(max(abs(astft)',minval)));
    colorbar;
    axis('xy');
    ylabel('Freq\div{}fs');
    v_cblabel('Mag (dB)');
    colormap(gca,rgb1);
    title(sprintf('Par:Off=%s,Scl=%s,Pad=%s,Del=%s', q.offset, q.scale, q.pad, q.groupdelay));
end
