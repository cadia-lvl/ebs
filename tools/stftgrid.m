function [stftg,metag]=stftgrid(stfte,meta,par)
% Interpolation of variable frame length STFT onto a regular grid
%
% Usage:    [stft,meta,grpd]=stfte(s,metain,[],par);                    % epoch-based STFT
%           [stft,meta]=stftgrid(stft,meta,par);                        % map onto a fixed grid unless par.interpstft='none'
%
%  Inputs: stfte(nframe,maxbin)     complex STFT coefficients
%          meta(nframe,nmeta)       metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples)]
%                                       although only the first three columns are used
%          maxfft                   max size of frame in stft domain [default: max(metain(:,2))]
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.interpstft  'none'      interpolation method in griddata: {'none','nearest','linear','natural','cubic','v4'}
%                                       par.interpHzps  1           Distance in frequency bins that is equivalent to a distance of one hop in time (bins per hop)
%                                       par.interpdom   'cplx'      Interpolate domain: {'cplx','magcph','crmcph'}
%                                       par.interpgrid  [...]       [nhop nbin] gives a time-frequency interpolation grid of nhop samples and fs/nbin Hz.
%                                                                   [nhop nbin] will be rounded to integers. If unspecified or empty, par.interpgrid = mean(meta(:,3)*[1 1].
%
% Outputs: stftg(nframe,nbin)       complex STFT coefficients
%          metag(nframe,nmeta)      output metadata with same column count as the input meta. metag(*,:)=[first-sample, frame-length=nbin, dft-length=nbin, 0, 1, 0]
%
% Algorithm Options:
%
%   par.interpstft  'none'          No STFT interpolation [default]
%                   'nearest'       Nearest-neighbour interpolation
%                   'linear'        Linear interpolation used in griddata.m
%                   'natural'       Uses Voronoi cells to determine neighbours for linear interpolation in griddata.m [1] [default]
%                   'cubic'         Cubic spline interpolation in griddata.m
%                   'v4'            MATLAB 4 method uses biharmonic spline interpolation in griddata.m [2]
%   par.interpdom   'cplx'          Interpolate in the complex domain
%                   'magcph'        Interpolate in magnitude domain and use phase from 'cplx' interpolation [default]
%                   'crmcph'        Interpolate in cube-root power domain and use phase from 'cplx' interpolation
%   par.interpgrid  [nhop nbin]     Time-frequency interpolation grid is nhop samples and fs/nbin Hz. Will be rounded to integers. [default = mean(metain(:,3)*[1 1]]
%   par.interpbph   bins_per_hop    Distance in frequency bins that is equivalent to a distance of one hop in time (bins per hop) [default: 1]
%                                       In different units, this is equivalent to bins_per_hop*fs^2/(nhop*nbin) Hz/second. A high value will tend
%                                       to smear the spectrogram in the frequency direction while a low value will smear it in the time direction.
%
% Refs: [1] Sibson, R. (1981). "A brief description of natural neighbor interpolation (Chapter 2)".
%           In V. Barnett. Interpreting Multivariate Data.  Chichester: John Wiley. pp. 21--36.
%       [2] David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data,
%           Geophysical Research Letters, 2, 139-142. 1987.
%
% Bugs/Suggestions:
% (1) Does not currently work well for very short frames (<7 samples)
% 
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.interpstft='natural';            % Interpolation method in griddata: {'none','nearest','linear','natural','cubic','v4'}
    q0.interpbph=1;                     % Distance in frequency bins that is equivalent to a distance of one hop time (bins per hop)
    q0.interpdom='magcph';              % Interpolatione domain: {'cplx','magcph','crmcph'}
    q0.interpgrid=[];                   % Time-frequency interpolation grid, [nhop nbin], is nhop samples and fs/nbin Hz. Will be rounded to integers. [default = mean(metain(:,3)*[1 1]]
end
%
% update algorithm parameters
if nargin<3
    q=q0;                                               % just copy default parameters
else
    q=v_paramsetch(q0,par);                             % update parameters from par input
end
if isempty(q.interpgrid)                                % no value specified for q.interpgrid
    q.interpgrid=repmat(round(mean(meta(:,3))),1,2);    % use the average DFT length
else
    q.interpgrid=round(q.interpgrid);                   % force integer values for the interpolation grid
end
if strcmp(q.interpstft,'none')
    stftg=stfte;
    metag=meta;
else
    [nframe,maxbin]=size(stfte);                        % number of frames and maximum fft size over all frames
    stftvv=stfte(:);                                    % make into a column vector
    taxi=meta(:,1:2)*[1;0.5]-0.5;                       % centre of input frames in samples (start @ 1)
    taxv=repmat(taxi,maxbin,1);                         % input centre-of-frame times in samples
    %%%% this doesn't work well for frames less than 7 samples long (luckily these are probably rare)
    nfftq=round(0.25*meta(:,3));                        % last quarter of of each row counts as negative frequencies
    faxv=reshape((mod(repmat(0:maxbin-1,nframe,1)+repmat(nfftq,1,maxbin),meta(:,3))-repmat(nfftq,1,maxbin))./repmat(meta(:,3),1,maxbin),[],1); % bins in fractions of fs
    % faxv=reshape((fs./meta(:,3))*(0:maxbin-1),[],1);
    msk=reshape(repmat(0:maxbin-1,nframe,1)<meta(:,3),[],1);    % first meta(:,3) entries in each row
    % msk=~isnan(stftvv);
    stftvv=stftvv(msk);                                 % eliminate non-existant entries from stft
    taxv=taxv(msk);                                     % ... and time coordinate in samples
    faxv=faxv(msk);                                     % and frequency coordinate in fractions of fs
    nhop=q.interpgrid(1);                               % frame hop in samples
    nbin=q.interpgrid(2);                               % effective fft length (not necessarily even)
    nbinp=1+floor(nbin/2);                              % number of positive output frequencies
    frst=max(ceil(taxi(1)-0.5*nbin+1.5),meta(1,1));     % start of first fixed frame in samples (frame-centre >= taxi(1)+1 and frame-start>=meta(1,1))
        % frst=max(ceil(taxi(1)-0.5*nhop+1.5),meta(1,1));     % start of first fixed frame in samples (frame-centre >=taxi(1)+1 to ensure no extrapolation at start)
    % frst=max(meta(1,1:2)*[1;0.5]-nbinp,1); % start of first fixed frame in samples (to ensure no extrapolation at start)
    nframef=floor(min((taxi(end)-frst-0.5-nbin*0.5)/nhop+1,(meta(end,1:2)*[1;1]-frst-nbin)/nhop+1));  % Num frames (frame-centre <= taxi(end)-1 and frame-end<=meta(end,1)+meta(end,2)-1)
    % nframef=1+floor(min((meta(end,1:2)*[2;1]-1-nbin-2*frst-1)/(2*nhop),(length(s)-frst-nbin+1)/nhop)); % #frames (to ensure no extrapolation at end)
    if nframef>0
        metarow=[nbin nbin 0 1 0];
        metag=[frst+(0:nframef-1)'*nhop repmat(metarow(1:min(max(size(meta,2),3),6)-1),nframef,1)];   % output metadata
        faxf=(0:nbinp-1)/nbin;                              % positive frequencies of fixed frame dft in fractions of fs
        taxf=metag(:,1)+nbin*0.5-0.5;            % centres of fixed frames in samples
        faxfv=reshape(repmat(faxf,nframef,1),[],1);         % fixed frame frequencies
        taxfv=repmat(taxf,nbinp,1);                         % fixed frame times
        taxfact=q.interpbph/(nbin*nhop);                    % relative weighting of time-frequency errors in fs/sample
        switch q.interpdom
            case 'cplx'                                     % interpolate complex values (i.e. real and imaginary separately)
                [xq,yq,vq]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                 % do complex interpolation onto fixed grid
            case 'magcph'                                   % interpolate magnitude and complex phase
                [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv),taxfv*taxfact,faxfv,q.interpstft);            % do magnitude interpolation onto fixed grid
                vqa=abs(vqc);                               % magnitude of complex interpolation
                msk=vqa~=0;                                 % complex phase irrelevant if vqa==0
                vq(msk)=vq(msk)./vqa(msk).*vqc(msk);        % magnitude from vq and phase from vqc unless vqc==0
            case 'crmcph'                                   % interpolate cube-root power and complex phase (as in Hermansky1990)
                [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv).^(2/3),taxfv*taxfact,faxfv,q.interpstft);     % do cube-root-power interpolation onto fixed grid
                vq=vq.^(3/2);                               % undo cube-root-power compression
                vqa=abs(vqc);                               % magnitude of complex interpolation
                msk=vqa~=0;                                 % complex phase irrelevant if vqa==0
                vq(msk)=vq(msk)./vqa(msk).*vqc(msk);        % magnitude from vq and phase from vqc unless vqc==0
        end
        stftg=reshape(vq,nframef,[]);                       % stft has one row per frame
        if mod(nbin,2)>0                                    % nbin is odd
            stftg=[real(stftg(:,1)) stftg(:,2:end) conj(stftg(:,end:-1:2))]; % force conjugate symmetry
        else                                                % nbin is even
            stftg=[real(stftg(:,1)) stftg(:,2:end-1) real(stftg(:,end)) conj(stftg(:,end-1:-1:2))]; % force conjugate symmetry
        end
    else % no frames to output
        stftg=zeros(0,nbin);
        metag=zeros(0,min(max(size(meta,2),3),6));
    end
end
if ~nargout && size(stftg,1)>0
    imagesc(taxf,faxf,db(abs(stftg(:,1:nbinp)))'); axis 'xy'; set(gca,'clim',get(gca,'clim')*[0;1]+[-40 0]);
    if ~strcmp(q.interpstft,'none')
        hold on
        plot(taxv,faxv,'.k');
        hold off
    end
    set(gca,'xlim',[meta(1,1) meta(end,1:2)*[1;1]-1],'ylim',[min(faxv)-0.02 max(faxv)+0.02]);
    colorbar;
    cblabel('dB');
    xlabel('Time (samples)');
    ylabel('Frequency / f_s');
end