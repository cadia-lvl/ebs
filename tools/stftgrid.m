function [stftg,metag]=stftgrid(stfte,meta,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine is now *OBSOLETE*, replace with:                                                                 %
%                                                                                                               %
%       par.interpfsps=par.interpbph/(par.interpgrid(1)*par.interpgrid(2));                                     %  
%       grid=[max(ceil(meta(1,1)+0.5*(meta(1,2)-1-par.interpgrid(2)+3)),meta(1,1)) par.interpgrid([2 1])];      %
%       [stft,meta]=stftegrid(stft,meta,grid,par);                                                              %
%                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of variable frame length STFT onto a regular grid
%
% Usage:    [stft,meta,grpd]=stfte(s,metain,[],par);                    % epoch-based STFT
%           [stft,meta]=stftgrid(stft,meta,par);                        % map onto a fixed grid unless par.interpstft='none'
%
%  Inputs: stfte(nframe,maxbin)     complex STFT coefficients
%          meta(nframe,nmeta)       metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples)]
%                                       although only columns 1,2,3 and 6 are currently used
%          maxfft                   max size of frame in stft domain [default: max(metain(:,2))]
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.interpstft  'none'      interpolation method in griddata: {'none','nearest','linear','natural','cubic','v4'}
%                                       par.interpbph   1           Distance in frequency bins that is equivalent to a distance of one hop in time (bins per hop)
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
%                   'natural'       Uses Voronoi cells to determine neighbours for linear interpolation in griddata.m [ref 1] [default]
%                   'cubic'         Cubic spline interpolation in griddata.m
%                   'v4'            MATLAB 4 method uses biharmonic spline interpolation in griddata.m [ref 2]
%   par.interpdom   'cplx'          Interpolate in the complex domain
%                   'magcph'        Interpolate in magnitude domain and use phase from 'cplx' interpolation [default]
%                   'crmcph'        Interpolate in cube-root power domain and use phase from 'cplx' interpolation
%   par.interpgrid  [nhop nbin]     Time-frequency interpolation grid is nhop samples and fs/nbin Hz. Will be rounded to integers. [default = mean(metain(:,3)*[1 1]]
%                                   can be calculated as [flen*fs/ovfact) 2*round(flen*fs/2)] where fs=sample freq (Hz), flen=equivalent frame length (s), ovfact=equivalent overlap factor
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
% (2) Code for calculating the group delay for each fixed frame does not take account of the data content
% (3) If fixed frame hop is large then some input frames may have no effect on the output (this may or may not be a problem)
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
if strcmp(q.interpstft,'none')                          % no interpolation requested ...
    stftg=stfte;                                        % ... copy across stft ...
    metag=meta;                                         % ... and meta data
else
    if size(meta,2)<6
        meta(1,6)=0;                                    % enlarge meta to 6 columns if necessary
    end
    [nframe,maxbin]=size(stfte);                        % number of frames and maximum fft size over all frames
    stftvv=stfte(:);                                    % make into a column vector
    taxi=meta(:,1:2)*[1;0.5]-0.5;                       % centre of input frames in samples (start @ 1)
    taxv=repmat(taxi,maxbin,1);                         % input centre-of-frame times in samples
    %%%% this doesn't work well for frames less than 7 samples long (luckily these are probably rare)
    nfftq=round(0.25*meta(:,3));                        % last quarter of of each row counts as negative frequencies
    faxv=reshape((mod(repmat(0:maxbin-1,nframe,1)+repmat(nfftq,1,maxbin),meta(:,3))-repmat(nfftq,1,maxbin))./repmat(meta(:,3),1,maxbin),[],1); % bins in fractions of fs
    msk=reshape(repmat(0:maxbin-1,nframe,1)<meta(:,3),[],1);    % first meta(:,3) entries in each row
    stftvv=stftvv(msk);                                 % eliminate non-existant entries from stft
    taxv=taxv(msk);                                     % ... and time coordinate in samples
    faxv=faxv(msk);                                     % ... and frequency coordinate in fractions of fs
    nhop=q.interpgrid(1);                               % frame hop in samples
    nbin=q.interpgrid(2);                               % effective fft length (not necessarily even)
    nbinp=1+floor(nbin/2);                              % number of positive output frequencies
    frst=max(ceil(taxi(1)-0.5*nbin+1.5),meta(1,1));     % start of first fixed frame in samples (frame-centre >= taxi(1)+1 and frame-start>=meta(1,1))
    nframef=floor(min((taxi(end)-frst-0.5-nbin*0.5)/nhop+1,(meta(end,1:2)*[1;1]-frst-nbin)/nhop+1));  % Num frames (frame-centre <= taxi(end)-1 and frame-end<=meta(end,1)+meta(end,2)-1)
    if nframef>0
        metag=[frst+(0:nframef-1)'*nhop repmat([nbin nbin 0 1 0],nframef,1)];   % output metadata
        faxf=(0:nbinp-1)/nbin;                              % positive frequencies of fixed frame dft in fractions of fs
        taxf=metag(:,1)+nbin*0.5-0.5;                       % centres of fixed frames in samples (start @ 1)
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
        % Calculate goup delay of output frames by using linear interpolation between the group delays of the input frames whose centres are either side of the
        % centre of the output frame while compensating for the starting sample of each of the frames.
        [xxi,xxf]=v_interval(taxf,taxi);                    % i'th fixed frame centre, taxf(i), lies between taxi(xxi(i)) and taxi(xxi(i)+1)
        msk=xxf>0.5;                                        % mask for fixed frames closer to xxi(i)+1 than to xxi(i)
        wtj=1-msk-(1-2*msk).*xxf;                           % weight to apply to gdfj below: 1-xxf(i) if msk(i)=0 or xxf(i) if mask(i)=1
        xxj=xxi(:)+msk;                                        % taxf(i) is closest to  taxi(xxj(i))
        xxk=xxi(:)+1-msk;                                      % other end of interval
        gdfj=v_modsym(meta(xxj,1)+meta(xxj,6),meta(xxj,3),taxf);    % add multple of DFT length to get assumed energy peak near centre of output frame 
        gdfk=v_modsym(meta(xxk,1)+meta(xxk,6),meta(xxk,3),gdfj);    % add multple of DFT length to get assumed energy peak near previous one 
        metag(:,6)=gdfj.*wtj+gdfk.*(1-wtj)-metag(:,1);      % group delay of fixed frame in samples: linear interpolate between gdfj and gdfk then compensate for start of frame
        if mod(nbin,2)>0                                    % nbin is odd
            stftg=[real(stftg(:,1)) stftg(:,2:end) conj(stftg(:,end:-1:2))]; % force conjugate symmetry
        else                                                % nbin is even
            stftg=[real(stftg(:,1)) stftg(:,2:end-1) real(stftg(:,end)) conj(stftg(:,end-1:-1:2))]; % force conjugate symmetry
        end
    else % no frames to output
        stftg=zeros(0,nbin);
        metag=zeros(0,6);
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