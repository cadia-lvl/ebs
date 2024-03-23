function [stft,meta]=stfte(s,metain,maxfft,par)
% Epoch-based stft
%
%  Inputs: s(ns,1)                  real-valued signal
%          metain(nframe,2)         input metadata: meta(*,:)=[first-sample,frame-length]
%          maxfft                   max size of frame in stft domain [default: max(metain(:,2))]
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.offset      'none'      offset removal: {'none','mean','ends'}
%                                       par.scale       'none'      scaling method: {'none','peakabs','rms'}
%                                       par.pad         'none'      zero-padding method: {'none','zero','ends'}
%                                       par.groupdelay  'none'      linear phase component: {'none','ewgd','cplx','phgr','gcif','gpdf','fmnb' + optional 'int' suffix}
%
% Outputs: stft(nframe,maxfft)      complex STFT coefficients
%          meta(nframe,6)           output metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay]
%
% Transformations are applied in the order offset, scale, pad, dft, groupdelay. If pad option is 'ends', the group delay can exceed the length of the unpadded frame.
%
% Algorithm Options:
%
%   par.offset      'none'      No offset-removal is performed so meta(:,4)=0
%                   'mean'      The mean of each frame is subtracted from the frame
%                   'ends'      The average of the first and last samples in each frame is subtracted from the frame
%   par.scale       'none'      No scaling is performed so meta(:,5)=1
%                   'peakabs'   Scale each frame so the peak absolute value equals 1
%                   'rms'       Scale each frame so that the root-mean-square value equals 1
%   par.pad         'none'      No padding is performed so dft-length for each frame equals the frame length and meta(:,3)=meta(:,2)
%                   'zero'      Pad each frame with zeros to a length of maxfft
%                   'ends'      Padd each frame with the average of the end values to a length of maxfft
%   par.groupdelay  'none'      No group delay compensation is performed so meta(:,6)=0 (i.e. zero samples delay)
%                   'ewgd'      Take group delay equal to the centre of gravity of the energy of the frame after padding
%                   'cplx'      As ewgd but convert position in frame to a complex number before doing the energy-weighted
%                               average. This mirrors the circularity of the DFT.
%                   'phgr'      Calculate the enegrgy-weighted phase decrement beween successive frequency bins in the DFT output
%                   'gcif'      Take group delay equal to par.gcifrac multiplied by the frame length, meta(:,2).
%                   'gpdf'      Take group delay equal to par.gpdfrac multiplied by the frame length, meta(:,2) [default=0.3].
%                   'fmnb'      Find optimum group delay subject to bounds par.fmbound as fraction of frame length [default bounds = [0.3 0.5]].
%                   '****int'   As above but rounded to an integer number of samples where '****' is one of the previous options
% Bugs/Suggestions:
% (1) for overlapping frames could apply a window before doing the fft
%
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.offset='none';                   %  offset removal: {'none','mean','ends'}
    q0.scale='none';                    %  scaling method: {'none','peakabs','rms'}
    q0.pad='none';                      % zero-padding method: {'none','zero','ends'}
    q0.groupdelay='none';               % linear phase component: {'none','ewgd','cplx','phgr','gcif','gpdf','fmnb' + optional 'int' suffix}
    q0.gpdfrac=0.3;                     % fraction of frame length for par.groupdelay='gpdf' option
    q0.fmbound=[0.3 0.5];               % bounds for group delay option par.groupdelay='fmbd' as fraction of frame length
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
meta=zeros(nframe,6);                   % space for metadata
meta(:,1:2)=metain;                     % copy frame start and frame length information
framelens=meta(:,2);                     % length of each frame in samples
if nargin<3 || isempty(maxfft)
    maxfft=max(framelens);              % set maxfft to maximum frame length
else
    framelens=min(framelens,maxfft);    % no framelength can exceed maxfft
end
sfr=zeros(nframe,maxfft);                                       % space for frames (one per row)
stft=NaN(nframe,maxfft);                                        % space for stft; unused entries will be left as NaN
sfrix=repmat(0:maxfft-1,nframe,1)+repmat(meta(:,1),1,maxfft);   % index into s for frame samples
sfrmk=repmat(1:maxfft,nframe,1)<=repmat(meta(:,2),1,maxfft);    % mask for valid values
sfr(sfrmk)=s(sfrix(sfrmk));                                     % enframe the data
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
    meta(:,5)=sqrt(sum(sfr.^2,2)./meta(:,2));               % find rms
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
elseif strcmp(q.scale,'peakabs')
    meta(:,5)=max(abs(sfr),[],2);                           % find max of abs(sfr)) in each frame
    sfr=sfr./repmat(meta(:,5)+(meta(:,5)==0),1,maxfft);     % scale the data (except when scale factor is zero)
else
    meta(:,5)=1;                                            % default scale factor is unity
end
if strcmp(q.pad,'none')                                     % if no padding
    meta(:,3)=meta(:,2);                                    % dft length equals the frame length
else                                                        % need to include padding
    meta(:,3)=maxfft;                                       % dft length equals maxfft
end
if isfield(par,'gcifrac')
    gcif=par.gcifrac;
else
    gcif=0;
end
if all(framelens==framelens(1))                                % all frames are the same length so we can do them all at once
    framelen=framelens(1);                                     % constant frame length for all frames
    if strcmp(q.pad,'ends') && framelen~=maxfft              % we need to pad frames with the average of the endpoints
        sfr(:,framelen+1:maxfft)=repmat(sfr(:,[1 framelen])*[0.5;0.5],1,maxfft-framelen);   % pad with average of endpoints
    end
    nfft=meta(1,3);                                         % DFT length of all frames
    stft(:,1:nfft)=fft(sfr(:,1:nfft),nfft,2);               % Perform DFT on all frames
    if ~strcmp(q.groupdelay,'none')
        switch q.groupdelay(1:4)
            case 'ewgd'
                meta(:,6)=sfr(:,1:nfft).^2*(0:nfft-1)'./sum(sfr(:,1:nfft).^2,2); % calculate EWGD for all frames
            case 'cplx'
                meta(:,6)=mod(angle(sfr(:,1:nfft).^2*exp(2i*pi*(0:nfft-1)'/nfft))*nfft/(2*pi),nfft); % calculate complex group delay for this frame
            case 'phgr'
                meta(:,6)=angle(sum(stft(:,1:nfft-1).*conj(stft(:,2:nfft)),2)./sum(abs(stft(:,1:nfft-1).*stft(:,2:nfft)),2))*nfft/(2*pi);
            case 'gcif'
                meta(:,6)=meta(:,2)*gcif;
            case 'gpdf'
                meta(:,6)=meta(:,2)*par.gpdfrac;
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
else                                                        % we must process frames individually 'cos varying lengths
    framelen=-1;                                            % initialize to invalid frame length (used to avoid recalculating window unnecessarily)
    for i=1:nframe
        framelen=framelens(i);                                 % length of this frame in samples
        if strcmp(q.pad,'ends') && framelen~=maxfft
            sfr(i,framelen+1:maxfft)=repmat(sfr(i,[1 framelen])*[0.5;0.5],1,maxfft-framelen); % pad with average of endpoints
        end
        nfft=meta(i,3);                                     % DFT length
        stft(i,1:nfft)=fft(sfr(i,1:nfft));                  % Perform DFT
        if ~strcmp(q.groupdelay,'none')
            switch q.groupdelay(1:4)
                case 'ewgd'
                    meta(i,6)=sfr(i,1:nfft).^2*(0:nfft-1)'/sum(sfr(i,1:nfft).^2); % calculate EWGD for this frame
                case 'cplx'
                    meta(i,6)=mod(angle(sfr(i,1:nfft).^2*exp(2i*pi*(0:nfft-1)'/nfft))*nfft/(2*pi),nfft); % calculate complex group delay for this frame
                case 'phgr'
                    meta(i,6)=mod(angle(stft(i,2:nfft-1)*stft(i,3:nfft)'/sum(abs(stft(i,2:nfft-1).*stft(i,3:nfft))))*nfft/(2*pi),nfft); % omit DC from calculation
                case 'gcif'
                    meta(i,6)=meta(i,2)*gcif;
                case 'gpdf'
                    meta(i,6)=meta(i,2)*par.gpdfrac;
                case 'fmnb'
                    % meta(i,6)=fminbnd(@(g) gderr(g,stft(i,1:nfft)),q.fmbound(1)*(framelen-1),q.fmbound(2)*(framelen-1)); % find optimum
                    [xx,yy,bb]=fftgdopterr(stft(i,1:nfft));                             % initial call computes fixed quantities for optimization
                    meta(i,6)= fminbnd(@(g) fftgdopterr(g,xx,yy,bb),q.fmbound(1)*(framelen-1),q.fmbound(2)*(framelen-1)); % find optimum
                otherwise
                    error(sprintf('par.goupdelay equals unknown value: %s',q.groupdelay));
            end
            if length(q.groupdelay)>4
                meta(i,6)=round(meta(i,6));                 % round group delay to an integer to an integer
            end
            stft(i,1:nfft)=stft(i,1:nfft).*exp(2i*pi/nfft*meta(i,6)*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1]); % apply group delay (except to Nyquist frequency)
        end
    end
end
if ~nargout
    rgb1=v_colormap('v_thermliny','k');                     % get colormap for magnitude
    rgb2=v_colormap('v_circrby','k');                       % get colormap for phase
    nfftp=1+floor(maxfft/2);                                % number of positive frequencies
    faxz=(0:nfftp-1)/maxfft;                                % frequency axis for spectrogram
    astft=zeros(nframe,nfftp);                      % space for the interpolated stft
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
