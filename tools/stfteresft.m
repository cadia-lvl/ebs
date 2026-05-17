function [stftg,metag]=stfteresft(stftx,meta,par,nfft,flen)
% resample epoch-based STFT in time and frequency
%
%  Inputs:  stftx(nfin,maxfft)  STFT coefficient array
%           meta(nfin,nmeta)    metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples)]
%           par                 parameters
%           nfft                output dft length
%           flen                output frame length in samples
%
%  Output:  stftg(nfout,nfft)   output STFT coefficient array
%           metag(nfout,nmeta)  output metadata with same column count as the input meta. metag(*,:)=[first-sample, frame-length=nbin, dft-length=nbin, 0, 1, 0]
%
%
%
persistent q0
%
% define default parameters
%
if isempty(q0)
q0.interpdom=  'magcph';       % Interpolatione domain: {'cplx','magcph','crmcph'}
q0.interpseq=  1;              % 0=interpolate T and F jointly; 1=interpolate T and F sequentially with an intermediate complex spectrum
q0.interpp=    5;              % length of intermediate filter
end
%
% update algorithm parameters
if nargin<3
    q=q0;                               % just copy default parameters
else
    q=v_paramsetch(q0,par);             % update parameters from par input
end
nfin=size(stftx,1);                         % number of inoput frames
nfout=floor((sum(meta(end,1:2))-meta(1,1))/flen); % number of output frames
stftc=zeros(nfin,nfft);                        % space for frequency-interpolated complex STFT
metag=[(0:nfout-1)'*flen+meta(1,1) repmat([flen nfft 1 0 0 0 0],nfout,1)];    % create output metadata with frame starting sample, frame len and DFT len
mtt=stftresm(meta(:,2),metag(:,2),q.interpp);
switch q.interpdom
    case 'cplx'                                         % interpolate complex values (i.e. real and imaginary separately)
        for i=1:nfin                                    % loop through each frame
            nfftx=meta(i,3);                            % input DFT length
            stftc(i,:)=stftx(i,1:nfftx)*fftresm(nfftx,nfft,q.interpp);   % do frequency interpolation
        end
        stftg=mtt*stftc;                                % do time interpolation
    case 'magcph'                                       % interpolate magnitude and complex phase
        stfty=zeros(nfin,nfft);                        % space for magnitude interpolated STFT
        for i=1:nfin                                    % loop through each frame
            nfftx=meta(i,3);                            % input DFT length
            mtf=fftresm(nfftx,nfft,q.interpp);       % frequency transformation matrix
            stftc(i,:)=stftx(i,1:nfftx)*mtf;            % interpolate complex STFT
            stfty(i,:)=abs(stftx(i,1:nfftx))*mtf;       % interpolate magnitude STFT
        end
        if q.interpseq                                % convert to intermediate complex STFT
            msk=stftc~=0;                               % only bother with phase if magnitude is non-zero
            stfty(msk)=stfty(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation (except where complex STFT is zero)
            stftc=mtt*stfty; % time-interpolate the complex intermediate STFT
            stftg=mtt*abs(stfty); % time-interpolate the magnitude
            msk=stftc~=0;  % mask for phase imposition
            stftg(msk)=stftg(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation (except where complex STFT is zero)
        else % leave intermediate STFT as magnitude+phase
            stftc=mtt*stftc; % time-interpolate the complex intermediate STFT
            stftg=mtt*stfty; % time-interpolate the magnitude
            msk=stftc~=0;  % mask for phase imposition
            stftg(msk)=stftg(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation (except where complex STFT is zero)
        end
    case 'crmcph'                                               % interpolate cube-root power and complex phase (as in Hermansky1990)
        stfty=zeros(nfin,nfft);                                % space for cube-root power frequency-interpolated STFT
        for i=1:nfin                                            % loop through each frame
            nfftx=meta(i,3);                                    % input DFT length
            mtf=fftresm(nfftx,nfft,q.interpp);               % frequency transformation matrix
            stftc(i,:)=stftx(i,1:nfftx)*mtf;                    % interpolate complex STFT
            stfty(i,:)=abs(stftx(i,1:nfftx)).^(2/3)*mtf;        % interpolate cube-root power STFT
        end
        if q.interpseq                                        % convert to intermediate complex STFT (with cube root power)
            msk=stftc~=0;                                       % only bother with phase if magnitude is non-zero
            stfty(msk)=stfty(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation on cube-root power
            stftc=mtt*stfty;                                    % time-interpolate the complex intermediate STFT
            stftg=(mtt*abs(stfty)).^(3/2);                      % time-interpolate the cube-root power and convert back to magnitude
            msk=stftc~=0;                                       % mask for phase imposition
            stftg(msk)=stftg(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation (except where complex STFT is zero)
        else                                                    % leave intermediate STFT as magnitude+phase
            stftc=mtt*stftc;                                    % time-interpolate the complex intermediate STFT
            stftg=(mtt*abs(stfty)).^(3/2);                      % time-interpolate the cube-root power and convert back to magnitude
            msk=stftc~=0;                                       % mask for phase imposition
            stftg(msk)=stftg(msk).*stftc(msk)./abs(stftc(msk)); % impose phase from complex interpolation (except where complex STFT is zero)
        end
end
