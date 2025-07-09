function [stftg,metag]=stftegrid(stfte,meta,grid,par)
% Interpolation of variable frame length STFT onto a regular grid
%
% Usage:    [stft,meta,grpd]=stfte(s,metain,[],par);                    % epoch-based STFT
%           [stft,meta]=stftegrid(stft,meta,grid,par);                        % map onto a fixed grid unless par.interpstft='none'
%
%  Inputs: stfte(nfin,maxbin)     complex STFT coefficients
%          meta(nfin,nmeta)       metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples)]
%                                       although only columns 1,2,3 and 6 are currently used
%          grid                     specifies the target grid in samples (starting @ 1): either [firstsamp(nfout) ndft(nfout)] for each frame or [firstsamp ndft nhop] for a uniform grid
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.interpstft  'none'      interpolation method in griddata: {'none','indep','nearest','linear','natural','cubic','v4'}
%                                       par.interpfsps  1e-5        fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second (fs per sample)
%                                       par.interpdom   'cplx'      Interpolation domain: {'cplx','magcph','crmcph'}
%                                       par.interpext   'rep'      Handling of extrapolated frames: {'omit','zero','rep','refl'}
%
% Outputs: stftg(nfout,nbin)       complex STFT coefficients
%          metag(nfout,nmeta)      output metadata with same column count as the input meta. metag(*,:)=[first-sample, frame-length=nbin, dft-length=nbin, 0, 1, 0]
%
% Algorithm Options:
%
%   par.interpstft  'none'          No STFT interpolation [default]
%                   'indep'         Interpolate in time and freq independently
%                   'nearest'       Nearest-neighbour interpolation
%                   'linear'        Linear interpolation used in griddata.m
%                   'natural'       Uses Voronoi cells to determine neighbours for linear interpolation in griddata.m [ref 1] [default]
%                   'cubic'         Cubic spline interpolation in griddata.m
%                   'v4'            MATLAB 4 method uses biharmonic spline interpolation in griddata.m [ref 2]
%   par.interpdom   'cplx'          Interpolate in the complex domain
%                   'magcph'        Interpolate in magnitude domain and use phase from 'cplx' interpolation [default]
%                   'crmcph'        Interpolate in cube-root power domain and use phase from 'cplx' interpolation
%   par.interfsps   fs/sample       fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second (fs per sample) [default: 1e-5]
%                                       A high value will tend to smear the spectrogram in the frequency direction while a low value will smear it in the time direction.
%   par.interpext   'omit'          Do not include any extrapolated frames [default]
%                   'zero'          Assume preceding/following input frames are zero
%                   'rep'           Assume preceding/following input frames replicate existing end frames
%                   'refl'          Reflect the preceding/following input frames
%
% Refs: [1] Sibson, R. (1981). "A brief description of natural neighbor interpolation (Chapter 2)".
%           In V. Barnett. Interpreting Multivariate Data.  Chichester: John Wiley. pp. 21--36.
%       [2] David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data,
%           Geophysical Research Letters, 2, 139-142. 1987.
%
% Bugs/Suggestions:
% (1) Does not currently work well for very short frames (<7 samples)
% (2) Code for calculating the group delay for each fixed frame does not take account of the data content or the length of the frame
% (3) If fixed frame hop is large then some input frames may have no effect on the output (this may or may not be a problem)
% (4) does not currently deal with scale factor and offset metadata; should apply these before interpolation and recalculate after
% (5) should have an option for doing time and frequency interpolation independently
% (6) if par.interpext='refl', should we include a replicate of the first/last frame + correctly replicate group delay information
% (7) when doint 2D interpolation, the calcuation of ninpre and ninpost should include any frames that might be in range for natural interpolation
% (8) might be more efficient to merge the three 1-D interpolation steps into a single stage
% (9) should perhaps force conjugate symmetry explicitly in 1-D case (e.g. in case DC and Nyquist are not exactly real)
%
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.interpstft='natural';            % Interpolation method in griddata: {'none','nearest','linear','natural','cubic','v4'}
    q0.interpfsps=1e-5;                     % Distance in frequency bins that is equivalent to a distance of one hop time (bins per hop)
    q0.interpdom='magcph';              % Interpolatione domain: {'cplx','magcph','crmcph'}
    q0.interpext=  'rep';              % handling of extrapolated frames
end
%
% update algorithm parameters
if nargin<3
    q=q0;                                               % just copy default parameters
else
    q=v_paramsetch(q0,par);                             % update parameters from par input
end
if strcmp(q.interpstft,'none')                              % no interpolation needed, so ...
    stftg=stfte;                                            % ... copy across stft ...
    metag=meta;                                             % ... and metadata
else                                                        % we need interpolation
    if size(meta,2)<6                                       % check if meta has all 6 columns
        if size(meta,2)<5
            meta(1,6)=0;                                    % enlarge meta to 6 columns
            meta(:,5)=1;                                    % and set gain to 1
        else
            meta(1,6)=0;                                    % enlarge meta to 6 columns
        end
    end
    [nfin,maxbin]=size(stfte);                              % number of input frames and maximum fft size over all frames
    finfix=all(meta(:,3)==meta(1,3));                       % true if input DFT length is fixed
    foutfix=all(grid(:,2)==grid(1,2));                      % true if output DFT length is fixed
    % sort out output grid
    if length(grid)==3                                      % grid=[firstsamp ndft nhop]
        nhop=grid(3);                                       % frame hop in samples
        nbin=grid(2);                                       % effective fft length (not necessarily even)
        frst=grid(1);
        nfout=floor(min((meta(end,1:2)*[1;0.5]-1-frst-nbin*0.5)/nhop+1,(meta(end,1:2)*[1;1]-frst-nbin)/nhop+1));  % Num frames (frame-centre <= taxin(end)-1 and frame-end<=meta(end,1)+meta(end,2)-1)
        metag=[frst+(0:nfout-1)'*nhop repmat([nbin nbin 0 1 0],nfout,1)];   % output metadata
        maxbinout=nbin;
    else                                                    % grid=[firstsamp(nfout) ndft(nfout)]
        nfout=size(grid,1);                                 % number of output frames
        metag=[grid(:,[1 2 2]) repmat([0 1 0],nfout,1)];    % output metadata
        maxbinout=max(grid(:,2));
    end
    stftg=NaN(0,maxbinout);                                 % zero-frame output STFT in case nfout is or becomes zero
    if nfout>0                                              % if there are any frames to output ...
        maxbinout=size(stftg,2);                            %   max number of DFT bins in output
        taxin=meta(:,1:2)*[1;0.5]-0.5;                      %   centre of input frames in samples (start @ 1)
        taxout=metag(:,1:2)*[1;0.5]-0.5;                     %   centres of output frames in samples (start @ 1)
        interpindep=strcmp(q.interpstft,'indep'); % doing 1D interpolation
        if strcmp(par.interpext,'omit')                % we may need to delete some frames
            msk=taxout<=taxin(1)-interpindep | taxout>=taxin(end)+interpindep;        % output frames that require extrapolation (less stringent for 1D)
            taxout(msk)=[];                                 % delete frame frame-centre list
            metag(msk,:)=[];                                % and from output metadata
            nfout=length(taxout);                           % revised number of output frames
        end
        if nfout>0                                          % check again if there are any frames to output
            stftg=NaN(nfout,maxbinout);                     % space for output STFT
            if interpindep                 % if independent interpolation in time and frequency ...
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1D Interpolation is performed in three stages:                            %
                %   (1) frequency interpolation onto an intermediate frequency grid (nbinx) %
                %   (2) time interpolation onto the desired time grid (taxout)              %
                %   (3) frequency interpolation onto the desired frequency grid (nfout)     %
                % Stages (1) and/or (3) may be omitted if they are redundant                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % convert input stft to the intepolation domain
                switch q.interpdom
                    case 'magcph'                               % interpolate magnitude and complex phase
                        stfty=abs(stfte);
                    case 'crmcph'                               % interpolate cube-root power and complex phase (as in Hermansky1990)
                        stfty=abs(stfte).^(2/3);
                    otherwise
                        stfty=[];
                end
                % Stage (1): interpolate in frequency onto intermediate grid
                switch 2*foutfix+finfix                         % switch according to input and/or output being fixed frame
                    case 0
                        nbinx=max(maxbin,size(stftg,2));        % input and output frame sizes both variable
                    case 1
                        nbinx=meta(1,3);                        % fixed input frame size;  output frame sizes variable
                    case 2
                        nbinx=grid(1,2);                        % variable input frame size;  fixed output frame size
                    case 3
                        nbinx=min(meta(1,3),grid(1,2));         % input and output frame sizes both fixed; use input or output resolution, whichever is smaller
                end
                if finfix && nbinx==meta(1,3)                   % input frequency resolution is unchanged
                    stftx=stfte(:,1:nbinx);                     % copy existing complex stft
                    if ~isempty(stfty)
                        stfty=stfty(:,1:nbinx);                 % and magnitude-domain version if it exists
                    end
                else % need to change the frequency resolution to nbinx bins
                    kkf=repmat((0:nbinx-1)/nbinx,nfin,1).*repmat(meta(:,3),1,nbinx); % fractional index into input frequency bins
                    kklo=floor(kkf);
                    kkf=kkf-kklo;
                    kkhi=mod(kklo+1,meta(:,3)); % make references to fs wrap around to 0
                    stftx=stfte(repmat((1:nfin)',1,nbinx)+nfin*kklo).*(1-kkf)+stfte(repmat((1:nfin)',1,nbinx)+nfin*kkhi).*kkf; % interpolate in frequency direction
                    if ~isempty(stfty)
                        stfty=stfty(repmat((1:nfin)',1,nbinx)+nfin*kklo).*(1-kkf)+stfty(repmat((1:nfin)',1,nbinx)+nfin*kkhi).*kkf; % interpolate in frequency direction
                    end
                end
                % Stage (2): interpolate in time onto target time grid (taxout)
                if taxout(1)<taxin(1) || taxout(end)>taxin(end)                 % we need to cope with extrapolation
                    ninpre=sum(2*taxin(1)-taxout(1)>taxin); % number of reflected inputs we need to pre-append
                    ninpost=sum(2*taxin(end)-taxout(end)<taxin); % number of reflected inputs we need to append
                    if ninpre==nfin || ninpost==nfin % we could extrapolate further
                        error('The number of extrapolated output frames cannot exceed the number of input frames');
                    end
                    taxin=[2*taxin(1)-taxin(ninpre+1:-1:2); taxin; 2*taxin(end)-taxin(end-1:-1:end-ninpost)]; % add tops and tails
                    [tti,ttf]=v_interval(taxout,taxin,'cC');
                    ttg=1-ttf; % weight for low end of interval
                    switch par.interpext
                        case 'zero'     %         Assume preceding/following input frames are zero
                            ttg(tti<=ninpre | tti>ninpre+nfin)=0; % treat extrapolated input points as zeros
                            ttf(tti<ninpre | tti>=ninpre+nfin)=0; % treat extrapolated input points as zeros
                            ttj=min(max(tti-ninpre+1,1),nfin);
                            tti=min(max(tti-ninpre,1),nfin);
                                     meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'rep'     %          Assume preceding/following input frames replicate existing first/last frames
                            ttj=min(max(tti-ninpre+1,1),nfin);
                            tti=min(max(tti-ninpre,1),nfin);
                                      meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'refl'     %    Reflect the preceding/following input frames ????????? should we include a replicate of the first/last frame ???????????
                            ttj=nfin-abs(1+abs(tti-ninpre)-nfin);
                            tti=nfin-abs(1+abs(tti-ninpre-1)-nfin);
                                 meta=[meta(ninpre+1:-1:2,:); meta; meta(end-1:-1:end-ninpost,:)];
                    end
                else % no extrapolation needed
                    [tti,ttf]=v_interval(taxout,taxin,'cC');
                    ttj=tti+1; % index for high end of interval
                    ttg=1-ttf; % weight for low end of interval
                end
                stftx=stftx(tti,:).*repmat(ttg,1,nbinx)+stftx(ttj,:).*repmat(ttf,1,nbinx); % perform complex-valued linear interpolation
                if ~isempty(stfty)
                    stfty=stfty(tti,:).*repmat(ttg,1,nbinx)+stfty(ttj,:).*repmat(ttf,1,nbinx); % perform complex-valued linear interpolation
                end

                % Stage (3): interpolate in frequency onto target frequency grid (possibly different for each frame)
                kkf=repmat((0:maxbinout-1)*nbinx,nfout,1)./repmat(metag(:,3),1,maxbinout);  % fractional index into intermediate frequency bins
                kklo=floor(kkf);                                                            % note this might exceed nbinx-1 for some frames
                kkf=kkf-kklo;                                                               % fractional part
                kkhi=mod(kklo+1,nbinx);                                                % make references to fs wrap around to 0
                stftx=stftx(repmat((1:nfout)',1,maxbinout)+nfout*mod(kklo,nbinx)).*(1-kkf)+stftx(repmat((1:nfout)',1,maxbinout)+nfout*kkhi).*kkf; % interpolate in frequency direction
                if ~isempty(stfty)
                    stftg=stfty(repmat((1:nfout)',1,maxbinout)+nfout*mod(kklo,nbinx)).*(1-kkf)+stfty(repmat((1:nfout)',1,maxbinout)+nfout*kkhi).*kkf; % interpolate in frequency direction
                end
                switch q.interpdom
                    case 'cplx'                                     % interpolate complex values (i.e. real and imaginary separately)
                        stftg=stftx; % perform complex-valued linear interpolation
                    case 'magcph'                                   % interpolate magnitude and complex phase
                        stfta=abs(stftx);                               % magnitude of complex interpolation
                        msk=stfta~=0;                                 % complex phase irrelevant if vqa==0
                        stftg(msk)=stftg(msk)./stfta(msk).*stftx(msk);        % magnitude from vq and phase from vqc unless vqc==0
                    case 'crmcph'                                   % interpolate cube-root power and complex phase (as in Hermansky1990)
                        stftg=stftg.^(3/2);                               % undo cube-root-power compression
                        stfta=abs(stftx);                                % magnitude of complex interpolation
                        msk=stfta~=0;                                 % complex phase irrelevant if vqa==0
                        stftg(msk)=stftg(msk)./stfta(msk).*stftx(msk);        % magnitude from vq and phase from vqc unless vqc==0
                end
                stftg(repmat(1:maxbinout,nfout,1)>repmat(metag(:,3),1,maxbinout))=NaN;      % set invalid entries to NaN
            else                                                    % perform 2D interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %   2-D Interpolation   %
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % sort out input data coordinates
                ninpre=sum(2*taxin(1)-taxout(1)>=taxin); % number of reflected inputs we need to pre-append
                ninpost=sum(2*taxin(end)-taxout(end)<=taxin); % number of reflected inputs we need to append
                if ninpre==nfin || ninpost==nfin
                    error('Output frame times cannot be extrapolated from input frame times');
                end
                if ninpre+ninpost>0 % we need to add reflected inputs
                    taxin=[2*taxin(1)-taxin(ninpre+1:-1:2); taxin; 2*taxin(end)-taxin(end-1:-1:end-ninpost)]; % add tops and tails to frame times
                    switch par.interpext
                        case 'zero'     %         Assume preceding/following input frames are zero
                            stfte=[zeros(ninpre,maxbin); stfte; zeros(ninpost,maxbin)];
                            meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'rep'     %          Assume preceding/following input frames replicate existing end frames
                            stfte=[repmat(stfte(1,:),ninpre,1); stfte; repmat(stfte(end,:),ninpost,1)];
                            meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'refl'     %    Reflect the preceding/following input frames
                            stfte=[stfte(ninpre+1:-1:2,:); stfte; stfte(end-1:-1:end-ninpost,:)];
                            meta=[meta(ninpre+1:-1:2,:); meta; meta(end-1:-1:end-ninpost,:)];
                    end
                    nfin=size(stfte,1); % update number of input frames
                end
                stftvv=stfte(:);                                    % make into a column vector
                taxv=repmat(taxin,maxbin,1);                        % input centre-of-frame times in samples (start @ 1)
                %%%% the next few lines don't work well for frames less than 7 samples long (luckily these are probably rare) but means we don't have to add frequency samples
                % better would be to calculate how many wrap-around samples we need at each end to ensure that the vertices needed for interpolation are present
                nfftq=round(0.25*meta(:,3));                        % last quarter of of each row counts as negative frequencies
                faxv=reshape((mod(repmat(0:maxbin-1,nfin,1)+repmat(nfftq,1,maxbin),meta(:,3))-repmat(nfftq,1,maxbin))./repmat(meta(:,3),1,maxbin),[],1); % bins in fractions of fs
                msk=reshape(repmat(0:maxbin-1,nfin,1)<meta(:,3),[],1);    % first meta(:,3) entries in each row
                stftvv=stftvv(msk);                                 % eliminate non-existant entries from stft
                taxv=taxv(msk);                                     % ... and time coordinate in samples
                faxv=faxv(msk);                                     % ... and frequency coordinate in fractions of fs
                % sort out output data coordinates
                taxfv=repmat(taxout,maxbinout,1);                   % output frame times
                faxfv=reshape(repmat(0:maxbinout-1,nfout,1)./repmat(metag(:,3),1,maxbinout),[],1);                      % bins in fractions of fs (may exceed 1 for some frames)
                mskp=repmat(0:maxbinout-1,nfout,1)<=repmat(0.5*metag(:,3),1,maxbinout);                                 % mask for positive frequencies
                taxfv=taxfv(mskp);
                faxfv=faxfv(mskp);
                taxfact=q.interpfsps;                                                                                   % relative weighting of time-frequency errors in fs/sample
                switch q.interpdom
                    case 'cplx'                                                                                         % interpolate complex values (i.e. real and imaginary separately)
                        [xq,yq,vq]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                 % do complex interpolation onto fixed grid
                    case 'magcph'                                                                                       % interpolate magnitude and complex phase
                        [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                        [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv),taxfv*taxfact,faxfv,q.interpstft);            % do magnitude interpolation onto fixed grid
                        vqa=abs(vqc);                                                                                   % magnitude of complex interpolation
                        msk=vqa~=0;                                                                                     % complex phase irrelevant if vqa==0
                        vq(msk)=vq(msk)./vqa(msk).*vqc(msk);                                                            % magnitude from vq and phase from vqc unless vqc==0
                    case 'crmcph'                                                                                       % interpolate cube-root power and complex phase (as in Hermansky1990)
                        [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                        [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv).^(2/3),taxfv*taxfact,faxfv,q.interpstft);     % do cube-root-power interpolation onto fixed grid
                        vq=vq.^(3/2);                                                                                   % undo cube-root-power compression
                        vqa=abs(vqc);                                                                                   % magnitude of complex interpolation
                        msk=vqa~=0;                                                                                     % complex phase irrelevant if vqa==0
                        vq(msk)=vq(msk)./vqa(msk).*vqc(msk);                                                            % magnitude from vq and phase from vqc unless vqc==0
                end
                stftg(mskp)=vq;                                                                                         % insert intepolated values into stftg
                for i=1:nfout                                                                                           % loop through frames to force conjugate symmetry
                    nbin=metag(i,3);                                                                                    % DFT length for this frame
                    nbinp=1+floor(nbin/2);                                                                              % number of positive frequencies
                    if mod(nbin,2)>0                                                                                    % nbin is odd
                        stftg(i,:)=[real(stftg(i,1)) stftg(i,2:nbinp) conj(stftg(i,nbinp:-1:2)) NaN(1,maxbinout-nbin)]; % force conjugate symmetry (odd nbin)
                    else                                                                                                % nbin is even
                        stftg(i,:)=[real(stftg(i,1)) stftg(i,2:nbinp-1) real(stftg(i,nbinp)) conj(stftg(i,nbinp-1:-1:2)) NaN(1,maxbinout-nbin)]; % force conjugate symmetry (even nbin)
                    end
                end
            end  % end of 1D or 2D interpolation                                                   % if strcmp(q.interpstft,'indep') ... else ... end
            % Calculate goup delay of output frames by using linear interpolation between the group delays of the input frames whose centres are either side of the
            % centre of the output frame while compensating for the starting sample of each of the frames.
            [xxi,xxf]=v_interval(taxout,taxin);                    % i'th fixed frame centre, taxout(i), lies between taxin(xxi(i)) and taxin(xxi(i)+1)
            msk=xxf>0.5;                                        % mask for fixed frames closer to xxi(i)+1 than to xxi(i)
            wtj=1-msk-(1-2*msk).*xxf;                           % weight to apply to gdfj below: 1-xxf(i) if msk(i)=0 or xxf(i) if mask(i)=1
            xxj=xxi+msk;                                        % taxout(i) is closest to  taxin(xxj(i))
            xxk=xxi+1-msk;                                      % other end of interval
            gdfj=v_modsym(meta(xxj,1)+meta(xxj,6),meta(xxj,3),taxout);    % add multple of DFT length to get assumed energy peak near centre of output frame
            gdfk=v_modsym(meta(xxk,1)+meta(xxk,6),meta(xxk,3),gdfj);    % add multple of DFT length to get assumed energy peak near previous one
            metag(:,6)=gdfj.*wtj+gdfk.*(1-wtj)-metag(:,1);      % group delay of fixed frame in samples: linear interpolate between gdfj and gdfk then compensate for start of frame
        end                                                         % if nfout>0 ... end; check if any output frames after deleting invalid output frames
    end                                                             % if nfout>0 ... end; initial check if any output frames
end                                                                 % if strcmp(q.interpstft,'none') ... else ... end; check if any interpolation to be done
if ~nargout && size(stftg,1)>0
    imagesc(taxout,faxf,db(abs(stftg(:,1:nbinp)))'); axis 'xy'; set(gca,'clim',get(gca,'clim')*[0;1]+[-40 0]);
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