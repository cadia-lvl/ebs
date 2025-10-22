function [stftg,metag]=stftegrid(stftv,meta,grid,par,lsym)
% Interpolation of variable frame length STFT onto a regular grid + multiple layers
%
% Usage:    [stftv,metav]=stfte(s,metain,[],par);                   % epoch-based STFT
%           [stft,meta]=stftegrid(stftv,metav,grid,par);            % map onto a fixed grid unless par.interpstft='none'
%
%  Inputs: stftv(nfin,maxbin,nlay)     STFT coefficient array with nlay layers
%          meta(nfin,nmeta)       metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples)]
%                                       although only columns 1,2,3 and 6 are currently used
%          grid                     specifies the target grid in samples (starting @ 1): either [firstsamp(nfout) framelen(nfout) ndft(nfout)] for each frame or [firstsamp framelen ndft nhop] for a uniform grid
%                                   NOTE: for backward compatibility [firstsamp(nfout) ndft(nfout)] for each frame or [firstsamp ndft nhop] for a uniform grid is also accepted
%          par                      parameter structure containing optional parameters
%                                       =Parameter=     =Default=   =Description=
%                                       par.interpdom=   'magcph';       % Interpolatione domain: {'cplx','magcph','crmcph'}
%                                       par.interpext=   'rep';          % Handling of extrapolated frames: {'omit','zero','rep','refl'}
%                                       par.interpstft=  'natural';      % Interpolation method for call to griddata: {'none','indep','nearest','linear','natural','cubic','v4'}
%                                       par.interpfsps=  1e-5;           % Dimensionless interpolation scale factor: fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second
%                                       par.interpgd=    'none';         % Interpolation of frame group delay: {'none','lin','linrep'}
%                                       par.interpof=    'none';         % Interpolation of frame offset: {'none','lin'}
%                                       par.interpsc=    'none';         % Interpolation of frame scale factor: {'none','lin','log'}
%                                       par.interptz=    'none';         % Compensation of phase for the frame starting sample: {'none','origin'} 
%          lsym(nlay)               Symmetry of layers: -2 for STFT or {0,1}={aligned,shifted} frequencies plus {0,2,4}= {Complex Hermitian, Real Symmetric, Real antisymmetric} [default: [-2; zeros(nlay-1,1)]]
%
%
% Outputs: stftg(nfout,nbin,nlay)  STFT coefficient array
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
%   par.interfsps   fs/sample       Dimensionless factor: fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second [default: 1e-5]
%                                       A high value will tend to smear the spectrogram in the frequency direction while a low value will smear it in the time direction.
%   par.interpext   'omit'          Do not include any extrapolated frames [default]
%                   'zero'          Assume preceding/following input frames are zero
%                   'rep'           Assume preceding/following input frames replicate existing end frames
%                   'refl'          Reflect the preceding/following input frames
%   par.interpgd    'none'          Group delay in output metadata will equal zero for all frames
%                   'lin'           Group delay in output metadata will be linearly interpolated based on frame centres while compensating for frame starting positions
%                   'linrep'        as for option 'lin' but input values will be interpreted modulo the DFT length
%   par.interpof    'none'          Offset in output metadata will equal zero for all frames
%                   'lin'           Offset in output metadata will be linearly interpolated based on frame centres
%   par.interpsc    'none'          Scale in output metadata will equal zero for all frames
%                   'lin'           Scale in output metadata will be linearly interpolated based on frame centres
%                   'log'           Scale in output metadata will be linearly interpolated in the log domain based on frame centres
%   par.interptz    'none'          Do not compensate phase for the frame starting sample
%                   'origin'        Refer all phases to to sample 0 before interpolation and then to frame start after interpolation
%
% Refs: [1] Sibson, R. (1981). "A brief description of natural neighbor interpolation (Chapter 2)".
%           In V. Barnett. Interpreting Multivariate Data.  Chichester: John Wiley. pp. 21--36.
%       [2] David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data,
%           Geophysical Research Letters, 2, 139-142. 1987.
%
% Bugs/Suggestions:
% (1) Does not currently work well for very short frames (<7 samples). We should instead calculate precisely which additional freqency bins are needed (similar to ninpre/ninpost)
% (2) Interpolation weights do not take account of the data content or the length of the frame (just the location of the frame centre). This might matter for the metadata
% (3) If fixed frame hop is large then some input frames may have no effect on the output (this may or may not be a problem)
% (6) if par.interpext='refl', should we include a replicate of the first/last frame + correctly replicate group delay information?
% (7) when doing 2D interpolation, the calcuation of ninpre and ninpost should include any frames that might be in range for natural interpolation
% (8) might be more efficient to merge the three 1-D interpolation steps into a single stage; however the interpolation would then need 8 coefficients rather than 2
% (9) should perhaps force conjugate symmetry explicitly in 1-D case (e.g. in case DC and Nyquist are not exactly real)
% (10) when adding extra frame at start and end (see ninpre/ninpost) we need to check that the metadata is corectly calculated
% (11) would it be more efficient to do 1D interpolation only for positive frequencies and then impose conjugate symmetry at the end (as for 2D)
% (12) Group delay compensation in line 145 is not right for Nyquist frequency if delay is an odd number of samples (not obvious what the solution is)
% (13) Group delay compensation in line 516 causes phase discontinuities if the frequency resolution has be made finer. e.g. if the frequency resolution is doubled then alternate output frequencies in alternate frames will be multiplied by -1.
% (14) DC and Nyquist entries are not handled quite right when enforcing antisymmetry. E.g. pi or N is antsymmetric with itself modulo 2pi or
% modulo N.
persistent q0
%
% define default parameters
%
if isempty(q0)
    q0.interpdom=   'magcph';       % Interpolatione domain: {'cplx','magcph','crmcph'}
    q0.interpext=   'rep';          % Handling of extrapolated frames: {'omit','zero','rep','refl'}
    q0.interpstft=  'natural';      % Interpolation method for call to griddata: {'none','indep','nearest','linear','natural','cubic','v4'}
    q0.interpfsps=  1e-5;           % Dimensionless interpolation scale factor: fs^-2 multiplied by the distance in Hz that is equivalent to a distance of one second
    q0.interpgd=    'none';         % Interpolation of frame group delay: {'none','lin','linrep'}
    q0.interpof=    'none';         % Interpolation of frame offset: {'none','lin'}
    q0.interpsc=    'none';         % Interpolation of frame scale factor: {'none','lin','log'}
    q0.interptz=    'none';         % Compensation of phase for the frame starting sample: {'none','origin'}
end
%
% update algorithm parameters
if nargin<3
    q=q0;                                                   % just copy default parameters
else
    q=v_paramsetch(q0,par);                                 % use the par input to update the parameters in q0
end
interptzq=1-strcmp(q.interptz,'none');                      % flag to indicate that we must compensate phase for the frame starting sample
if strcmp(q.interpstft,'none')                              % no interpolation needed, so ...
    stftg=stftv;                                            % ... copy across stft ...
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
    [nfin,maxbin,nlay]=size(stftv);                         % number of input frames and maximum fft size over all frames
    if nargin<5
        lsym=[-2;zeros(nlay-1,1)];                                 % default is Complex Hermitian symmetry for all layers
    end
    finfix=all(meta(:,3)==meta(1,3));                       % true if input DFT length is fixed
    % sort out output grid
    [nfout,grcol]=size(grid);                                 % number of output frames
    grcol3=grcol==3; % flag for 3 columns
    if nfout==1 && (grcol==3 || grcol==4)                     % grid = [firstsamp framelen ndft nhop] or, legacy only, [firstsamp ndft nhop]
        nhop=grid(4-grcol3);                                       % frame hop in samples
        nbin=grid(3-grcol3);                                       % effective fft length (not necessarily even)
        frlen=grid(2);                                      % frame length
        frst=grid(1);                                       % start of first frame
        nfout=floor(min((meta(end,1:2)*[1;0.5]-1-frst-frlen*0.5)/nhop+1,(meta(end,1:2)*[1;1]-frst-frlen)/nhop+1));  % Num frames (frame-centre <= taxin(end)-1 and frame-end<=meta(end,1)+meta(end,2)-1)
        metag=[frst+(0:nfout-1)'*nhop repmat([frlen nbin 0 1 0],nfout,1)];   % output metadata
        maxbinout=nbin;
    else                                                    % grid = [firstsamp(nfout) framelen(nfout) ndft(nfout)] or, legacy only, [firstsamp(nfout) ndft(nfout)]
        nfout=size(grid,1);                                 % number of output frames
        metag=[grid(:,[1 2 2+grcol3]) repmat([0 1 0],nfout,1)];    % initialize output metadata
        maxbinout=max(grid(:,3));
    end
    foutfix=all(metag(:,3)==metag(1,3));                    % true if output DFT length is fixed
    stftg=NaN(0,maxbinout,nlay);                                 % zero-frame output STFT in case nfout is or becomes zero
    if nfout>0                                              % if there are any frames to output ...
        maxbinout=size(stftg,2);                            %   max number of DFT bins in output
        taxin=meta(:,1:2)*[1;0.5]-0.5;                      %   centre of input frames in samples (start @ 1)
        taxout=metag(:,1:2)*[1;0.5]-0.5;                    %   centres of output frames in samples (start @ 1)
        interpindep=strcmp(q.interpstft,'indep');           % doing 1D interpolation
        if strcmp(par.interpext,'omit')                     % we may need to delete some frames
            msk=taxout<=taxin(1)-interpindep | taxout>=taxin(end)+interpindep;        % output frames that require extrapolation (less stringent for 1D)
            taxout(msk)=[];                                 % delete frame frame-centre list
            metag(msk,:)=[];                                % and from output metadata
            nfout=length(taxout);                           % revised number of output frames
        end
        if nfout>0                                          % check again if there are any frames to output
            stftg=NaN(nfout,maxbinout,nlay);                     % space for output STFT
            if interpindep                 % if independent interpolation in time and frequency ...
                for ilay=1:nlay
                    stftvl=zeros(nfin,maxbin); % space for this layer
                    if lsym(ilay)<0 % complex STFT so adjust group delay using metadata
                        for i=1:nfin % for now undo the group delay frame by frame compensating for start-sample, scale, and group-delay: meta(:,[1 5 6]).
                            nfft=meta(i,3);                     % DFT length of this frame
                            stftvl(i,1:nfft)=stftv(i,1:nfft,ilay).*exp(-2i*pi/nfft*(interptzq*meta(i,1)+meta(i,6))*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1])*meta(i,5); % apply non-integer group delay (except to Nyquist frequency)
                            stftvl(i,1)=stftvl(i,1)+meta(i,4)*meta(i,3);   % add offset*DFT_length
                        end
                    else
                        stftvl=stftv(:,:,ilay); % extract the current layer
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % 1D Interpolation is performed in three stages:                            %
                    %   (1) frequency interpolation onto an intermediate frequency grid (nbinx) %
                    %   (2) time interpolation onto the desired time grid (taxout)              %
                    %   (3) frequency interpolation onto the desired frequency grid (nfout)     %
                    % Stages (1) and/or (3) may be omitted if they are redundant                %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if lsym(ilay)<2                                                                                             % complex hermitian data
                        % convert input stft to the intepolation domain
                        switch q.interpdom
                            case 'magcph'                               % interpolate magnitude and complex phase
                                stfty=abs(stftvl);
                            case 'crmcph'                               % interpolate cube-root power and complex phase (as in Hermansky1990)
                                stfty=abs(stftvl).^(2/3);
                            otherwise
                                stfty=[];
                        end
                        % Stage (1): interpolate in frequency onto intermediate grid
                        switch 2*foutfix+finfix                         % switch according to input and/or output being fixed frame
                            case 0
                                nbinx=max(maxbin,maxbinout);        % input and output frame sizes both variable: choose densest
                            case 1
                                nbinx=meta(1,3);                        % fixed input frame size;  output frame sizes variable
                            case 2
                                nbinx=metag(1,3);                       % variable input frame size;  fixed output frame size
                            case 3
                                nbinx=min(meta(1,3),metag(1,3));        % input and output frame sizes both fixed; use input or output resolution, whichever is smaller
                        end
                        if finfix && nbinx==meta(1,3)                   % input frequency resolution is unchanged
                            stftx=stftvl(:,1:nbinx);                     % copy existing complex stft
                            if ~isempty(stfty)
                                stfty=stfty(:,1:nbinx);                 % and magnitude-domain version if it exists
                            end
                        else                                                                                % need to change the frequency resolution to nbinx bins
                            if mod(lsym(ilay),2)                                                            % frequencies of this layer are shifted by 0.5 * frequency increment
                                kkf=repmat((0.5:nbinx-0.5)/nbinx,nfin,1).*repmat(meta(:,3),1,nbinx)-0.5;    % non-integer index into input frequency bins size=(nfin,nbinx) [base=0]
                            else                                                                            % frequencies of this layer are aligned with DFT
                                kkf=repmat((0:nbinx-1)/nbinx,nfin,1).*repmat(meta(:,3),1,nbinx);            % non-integer index into input frequency bins size=(nfin,nbinx) [base=0]
                            end
                            kklo=floor(kkf);
                            kkf=kkf-kklo;
                            kkhi=mod(kklo+1,meta(:,3));                 % make references to fs wrap around to 0
                            frix=repmat((1:nfin)',1,nbinx);             % portion of linear index due to frame
                            stftx=stftvl(frix+nfin*kklo).*(1-kkf)+stftvl(frix+nfin*kkhi).*kkf; % interpolate in frequency direction
                            if ~isempty(stfty)
                                stfty=stfty(frix+nfin*kklo).*(1-kkf)+stfty(frix+nfin*kkhi).*kkf; % interpolate in frequency direction
                            end
                        end
                        % Stage (2): interpolate in time onto target time grid (taxout)
                        if taxout(1)<taxin(1) || taxout(end)>taxin(end)     % we need to cope with extrapolation
                            ninpre=sum(2*taxin(1)-taxout(1)>taxin);         % number of reflected inputs we need to pre-append
                            ninpost=sum(2*taxin(end)-taxout(end)<taxin);    % number of reflected inputs we need to append
                            if ninpre==nfin || ninpost==nfin                % we could extrapolate further
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
                            stfty=stfty(tti,:).*repmat(ttg,1,nbinx)+stfty(ttj,:).*repmat(ttf,1,nbinx); % perform linear interpolation in chosen magnitude domain
                        end

                        % Stage (3): interpolate in frequency onto target frequency grid (possibly different for each frame)
                        if mod(lsym(ilay),2)                                                                                    % frequencies of this layer are shifted by 0.5 * frequency increment
                            kkf=repmat((0.5:maxbinout-0.5)*nbinx,nfout,1)./repmat(metag(:,3),1,maxbinout)-0.5;  % fractional index into intermediate frequency bins size=(nfout,maxbinout) [base=0]
                        else  % frequencies of this layer are aligned with DFT
                            kkf=repmat((0:maxbinout-1)*nbinx,nfout,1)./repmat(metag(:,3),1,maxbinout);  % fractional index into intermediate frequency bins size=(nfout,maxbinout) [base=0]
                        end
                        kklo=floor(kkf);                                                            % note this might exceed nbinx-1 for some frames
                        kkf=kkf-kklo;                                                               % fractional part
                        kkhi=mod(kklo+1,nbinx);                                                % make references to fs wrap around to 0
                        stftx=stftx(repmat((1:nfout)',1,maxbinout)+nfout*mod(kklo,nbinx)).*(1-kkf)+stftx(repmat((1:nfout)',1,maxbinout)+nfout*kkhi).*kkf; % interpolate complex stft in frequency direction
                        if ~isempty(stfty)
                            stftgl=stfty(repmat((1:nfout)',1,maxbinout)+nfout*mod(kklo,nbinx)).*(1-kkf)+stfty(repmat((1:nfout)',1,maxbinout)+nfout*kkhi).*kkf; % interpolate magnitude in frequency direction
                        end
                        switch q.interpdom
                            case 'cplx'                                         % interpolate complex values (i.e. real and imaginary separately)
                                stftgl=stftx;                                    % perform complex-valued linear interpolation
                            case 'magcph'                                       % interpolate magnitude and complex phase
                                stfta=abs(stftx);                               % magnitude of complex interpolation
                                msk=stfta~=0;                                   % complex phase irrelevant if vqa==0
                                stftgl(msk)=stftgl(msk)./stfta(msk).*stftx(msk);  % magnitude from vq and phase from vqc unless vqc==0
                            case 'crmcph'                                       % interpolate cube-root power and complex phase (as in Hermansky1990)
                                stftgl=stftgl.^(3/2);                             % undo cube-root-power compression
                                stfta=abs(stftx);                               % magnitude of complex interpolation
                                msk=stfta~=0;                                   % complex phase irrelevant if vqa==0
                                stftgl(msk)=stftgl(msk)./stfta(msk).*stftx(msk);  % magnitude from vq and phase from vqc unless vqc==0
                        end
                        stftgl(repmat(1:maxbinout,nfout,1)>repmat(metag(:,3),1,maxbinout))=NaN;      % set invalid entries to NaN
                        stftg(:,:,ilay)=stftgl; % insert into output array
                    else % real symmetric or antisymmetric data
                        % Stage (1): interpolate in frequency onto intermediate grid
                        switch 2*foutfix+finfix                         % switch according to input and/or output being fixed frame
                            case 0
                                nbinx=max(maxbin,maxbinout);        % input and output frame sizes both variable: choose densest
                            case 1
                                nbinx=meta(1,3);                        % fixed input frame size;  output frame sizes variable
                            case 2
                                nbinx=metag(1,3);                       % variable input frame size;  fixed output frame size
                            case 3
                                nbinx=min(meta(1,3),metag(1,3));        % input and output frame sizes both fixed; use input or output resolution, whichever is smaller
                        end
                        if finfix && nbinx==meta(1,3)                   % input frequency resolution is unchanged
                            stftx=stftvl(:,1:nbinx);                    % copy existing data
                        else                                            % need to change the frequency resolution to nbinx bins
                            if mod(lsym(ilay),2)                        % frequencies of this layer are shifted by 0.5 * frequency increment
                                kkf=mod(repmat((0.5:nbinx-0.5)/nbinx,nfin,1).*repmat(meta(:,3),1,nbinx)-0.5,meta(:,3));    % non-integer index into input frequency bins size=(nfin,nbinx) [base=0]
                            else                                        % frequencies of this layer are aligned with DFT
                                kkf=repmat((0:nbinx-1)/nbinx,nfin,1).*repmat(meta(:,3),1,nbinx);            % non-integer index into input frequency bins size=(nfin,nbinx) [base=0]
                            end
                            kklo=floor(kkf);
                            kkf=kkf-kklo;
                            kkhi=mod(kklo+1,meta(:,3));                 % make references to fs wrap around to 0
                            frix=repmat((1:nfin)',1,nbinx);             % portion of linear index due to frame
                            stftx=stftvl(frix+nfin*kklo).*(1-kkf)+stftvl(frix+nfin*kkhi).*kkf; % interpolate in frequency direction
                        end
                        % Stage (2): interpolate in time onto target time grid (taxout)
                        if taxout(1)<taxin(1) || taxout(end)>taxin(end)     % we need to cope with extrapolation
                            ninpre=sum(2*taxin(1)-taxout(1)>taxin);         % number of reflected inputs we need to pre-append
                            ninpost=sum(2*taxin(end)-taxout(end)<taxin);    % number of reflected inputs we need to append
                            if ninpre==nfin || ninpost==nfin                % we could extrapolate further
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
                        stftx=stftx(tti,:).*repmat(ttg,1,nbinx)+stftx(ttj,:).*repmat(ttf,1,nbinx); % perform linear interpolation

                        % Stage (3): interpolate in frequency onto target frequency grid (possibly different for each frame)
                        if mod(lsym(ilay),2)                                                                                    % frequencies of this layer are shifted by 0.5 * frequency increment
                            kkf=repmat((0.5:maxbinout-0.5)*nbinx,nfout,1)./repmat(metag(:,3),1,maxbinout)-0.5;  % fractional index into intermediate frequency bins size=(nfout,maxbinout) [base=0]
                        else  % frequencies of this layer are aligned with DFT
                            kkf=repmat((0:maxbinout-1)*nbinx,nfout,1)./repmat(metag(:,3),1,maxbinout);  % fractional index into intermediate frequency bins size=(nfout,maxbinout) [base=0]
                        end
                        kklo=floor(kkf);                                                            % note this might exceed nbinx-1 for some frames
                        kkf=kkf-kklo;                                                               % fractional part
                        kkhi=mod(kklo+1,nbinx);                                                % make references to fs wrap around to 0
                        stftgl=stftx(repmat((1:nfout)',1,maxbinout)+nfout*mod(kklo,nbinx)).*(1-kkf)+stftx(repmat((1:nfout)',1,maxbinout)+nfout*kkhi).*kkf; % interpolate complex stft in frequency direction
                        stftgl(repmat(1:maxbinout,nfout,1)>repmat(metag(:,3),1,maxbinout))=NaN;      % set invalid entries to NaN
                        stftg(:,:,ilay)=stftgl; % insert into output array
                    end
                end % loop for each layer
            else                                                        % perform 2D interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   2-D Interpolation                                                                        %
                %                                                                                            %
                %  Note that we only interpolate the first half of the spectrum (i.e. +ve frequencies) and   %
                %  use symmetry to calculate the remainder                                      %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sort out input data coordinates
                ninpre=sum(2*taxin(1)-taxout(1)>=taxin); % number of reflected inputs we need to pre-append
                ninpost=sum(2*taxin(end)-taxout(end)<=taxin); % number of reflected inputs we need to append
                if ninpre==nfin || ninpost==nfin
                    error('Output frame times cannot be extrapolated from input frame times');
                end
                if ninpre+ninpost>0 % we need to add reflected inputs in time dimension
                    taxin=[2*taxin(1)-taxin(ninpre+1:-1:2); taxin; 2*taxin(end)-taxin(end-1:-1:end-ninpost)]; % add tops and tails to frame times
                    switch par.interpext
                        case 'zero'                                 %         Assume preceding/following input frames are zero
                            stftv=cat(1,zeros(ninpre,maxbin,nlay),stftv,zeros(ninpost,maxbin,nlay));
                            meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'rep'                                  %          Assume preceding/following input frames replicate existing end frames
                            stftv=cat(1,repmat(stftv(1,:,:),[ninpre,1,1]),stftv,repmat(stftv(end,:,:),[ninpost,1,1]));
                            meta=[repmat(meta(1,:),ninpre,1); meta; repmat(meta(end,:),ninpost,1)];
                        case 'refl'                                 %    Reflect the preceding/following input frames
                            stftv=cat(1,stftv(ninpre+1:-1:2,:,:),stftv,stftv(end-1:-1:end-ninpost,:,:));
                            meta=[meta(ninpre+1:-1:2,:); meta; meta(end-1:-1:end-ninpost,:)];
                    end
                    nfin=size(stftv,1);                             % update number of input frames to include the added frames
                end
                taxv=repmat(taxin,maxbin,1);                        % input centre-of-frame times in samples (start @ 1)
                nfftq=round(0.25*meta(:,3));                        % last quarter of of each row counts as negative frequencies
                msk=reshape(repmat(0:maxbin-1,nfin,1)<meta(:,3),[],1);    % first meta(:,3) entries in each row (as a column vector mask)
                taxv=taxv(msk);                                     % eliminate non-existant entries time coordinate in samples
                mskp=repmat(0:maxbinout-1,nfout,1)<=repmat(0.5*metag(:,3),1,maxbinout);                                 % mask to restrict output to positive frequencies only
                taxfv=repmat(taxout,maxbinout,1);                   % output frame times
                taxfv=taxfv(mskp);
                taxfact=q.interpfsps;                               % relative weighting of time-frequency errors in fs/sample
                stftgl=NaN(nfout,maxbinout);                        % space for output STFT for this layer
                for ilay=1:nlay
                    stftvl=zeros(nfin,maxbin);                      % space for this layer
                    if lsym(ilay)<0                                 % complex STFT so adjust group delay using metadata
                        for i=1:nfin                                % for now undo the group delay frame by frame
                            nfft=meta(i,3);                         % DFT length of this frame
                            stftvl(i,1:nfft)=stftv(i,1:nfft,ilay).*exp(-2i*pi/nfft*(interptzq*meta(i,1)+meta(i,6))*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1])*meta(i,5); % apply non-integer group delay (except to Nyquist frequency)
                            stftvl(i,1)=stftvl(i,1)+meta(i,4)*meta(i,3);            % add offset*DFT_length
                        end
                        stftvv=stftvl(:);                                           % make current layer into a column vector
                    else
                        stftvv=reshape(stftv(:,:,ilay),[],1);                       % make current layer into a column vector
                    end
                    %%%% the next few lines don't work well for frames less than 7 samples long (luckily these are probably rare) but means we don't have to add frequency samples
                    % better would be to calculate how many wrap-around samples we need at each end to ensure that the vertices needed for interpolation are present
                    if mod(lsym(ilay),2)                                                                        % frequencies of this layer are shifted by 0.5 * frequency increment
                        faxv=reshape((mod(repmat(0.5:maxbin-0.5,nfin,1)+repmat(nfftq,1,maxbin),meta(:,3))-repmat(nfftq,1,maxbin))./repmat(meta(:,3),1,maxbin),[],1); % bins in fractions of fs
                        faxfv=reshape(repmat(0.5:maxbinout-0.5,nfout,1)./repmat(metag(:,3),1,maxbinout),[],1);  % bins in fractions of fs (may exceed 1 for some frames)
                    else                                                                                        % frequencies of this layer are aligned with STFT
                        faxv=reshape((mod(repmat(0:maxbin-1,nfin,1)+repmat(nfftq,1,maxbin),meta(:,3))-repmat(nfftq,1,maxbin))./repmat(meta(:,3),1,maxbin),[],1); % bins in fractions of fs
                        faxfv=reshape(repmat(0:maxbinout-1,nfout,1)./repmat(metag(:,3),1,maxbinout),[],1);      % bins in fractions of fs (may exceed 1 for some frames)
                    end
                    stftvv=stftvv(msk);                                                                         % eliminate non-existant entries from stft
                    faxv=faxv(msk);                                                                             % ... and frequency coordinate in fractions of fs
                    % sort out output data coordinates
                    faxfv=faxfv(mskp);                                                                          % restrict output interpolation positions to positive frequencies only
                    if lsym(ilay)<2                                                                                             % complex hermitian data
                        switch q.interpdom
                            case 'cplx'                                                                                         % interpolate complex values (i.e. real and imaginary separately)
                                [xq,yq,vq]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                 % do complex interpolation onto fixed grid
                            case 'magcph'                                                                                       % interpolate magnitude and complex phase
                                [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                                [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv),taxfv*taxfact,faxfv,q.interpstft);            % do magnitude interpolation onto fixed grid
                                vqa=abs(vqc);                                                                                   % magnitude of complex interpolation
                                mskv=vqa~=0;                                                                                     % complex phase irrelevant if vqa==0
                                vq(mskv)=vq(mskv)./vqa(mskv).*vqc(mskv);                                                            % magnitude from vq and phase from vqc unless vqc==0
                            case 'crmcph'                                                                                       % interpolate cube-root power and complex phase (as in Hermansky1990)
                                [xq,yq,vqc]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                % do complex interpolation onto fixed grid
                                [xq,yq,vq]=griddata(taxv*taxfact,faxv,abs(stftvv).^(2/3),taxfv*taxfact,faxfv,q.interpstft);     % do cube-root-power interpolation onto fixed grid
                                vq=vq.^(3/2);                                                                                   % undo cube-root-power compression
                                vqa=abs(vqc);                                                                                   % magnitude of complex interpolation
                                mskv=vqa~=0;                                                                                     % complex phase irrelevant if vqa==0
                                vq(mskv)=vq(mskv)./vqa(mskv).*vqc(mskv);                                                            % magnitude from vq and phase from vqc unless vqc==0
                        end
                        stftgl(mskp)=vq;                                                                                            % insert interpolated values into stftg layer
                        if mod(lsym(ilay),2)                                                                                        % frequencies of this layer are shifted by 0.5 * frequency increment
                            for i=1:nfout                                                                                           % loop through frames to force conjugate symmetry
                                nbin=metag(i,3);                                                                                    % DFT length for this frame
                                nbinp=1+floor(nbin/2);                                                                              % number of positive frequencies
                                if mod(nbin,2)>0                                                                                    % nbin is odd
                                    stftg(i,:,ilay)=[stftgl(i,1:nbinp-1) real(stftgl(i,nbinp)) conj(stftgl(i,nbinp-1:-1:1)) NaN(1,maxbinout-nbin)]; % force conjugate symmetry (odd nbin, shifted)
                                else                                                                                                % nbin is even
                                    stftg(i,:,ilay)=[stftgl(i,1:nbinp-1) conj(stftgl(i,nbinp-1:-1:1)) NaN(1,maxbinout-nbin)];       % force conjugate symmetry (even nbin, shifted)
                                end
                            end
                        else                                                                                                        % frequencies of this layer are aligned with STFT
                            for i=1:nfout                                                                                           % loop through frames to force conjugate symmetry
                                nbin=metag(i,3);                                                                                    % DFT length for this frame
                                nbinp=1+floor(nbin/2);                                                                              % number of positive frequencies
                                if mod(nbin,2)>0                                                                                    % nbin is odd
                                    stftg(i,:,ilay)=[real(stftgl(i,1)) stftgl(i,2:nbinp) conj(stftgl(i,nbinp:-1:2)) NaN(1,maxbinout-nbin)]; % force conjugate symmetry (odd nbin, unshifted)
                                else                                                                                                % nbin is even
                                    stftg(i,:,ilay)=[real(stftgl(i,1)) stftgl(i,2:nbinp-1) real(stftgl(i,nbinp)) conj(stftgl(i,nbinp-1:-1:2)) NaN(1,maxbinout-nbin)]; % force conjugate symmetry (even nbin, unshifted)
                                end
                            end
                        end
                    else                                                                                                            % lsym(ilay)>=2 so data is real symmetric or antisymmetric
                        [xq,yq,vq]=griddata(taxv*taxfact,faxv,stftvv,taxfv*taxfact,faxfv,q.interpstft);                             % do linear interpolation onto fixed grid
                        stftgl(mskp)=vq;                                                                                            % insert interpolated values into stftg layer
                        if lsym*ilay>=4                                                                                             % real antisymmnetric
                            if mod(lsym(ilay),2)                                                                                    % frequencies of this layer are shifted by 0.5 * frequency increment
                                for i=1:nfout                                                                                       % loop through frames to force conjugate symmetry
                                    nbin=metag(i,3);                                                                                % DFT length for this frame
                                    nbinp=1+floor(nbin/2);                                                                          % number of positive frequencies
                                    if mod(nbin,2)>0                                                                                % nbin is odd
                                        stftg(i,:,ilay)=[stftgl(i,1:nbinp-1) 0 -stftgl(i,nbinp-1:-1:1) NaN(1,maxbinout-nbin)];      % force antisymmetry (odd nbin)
                                    else                                                                                            % nbin is even
                                        stftg(i,:,ilay)=[stftgl(i,1:nbinp-1) -stftgl(i,nbinp-1:-1:1) NaN(1,maxbinout-nbin)];        % force antisymmetry (even nbin)
                                    end
                                end
                            else                                                                                                    % frequencies of this layer are aligned with STFT
                                for i=1:nfout                                                                                       % loop through frames to force conjugate symmetry
                                    nbin=metag(i,3);                                                                                % DFT length for this frame
                                    nbinp=1+floor(nbin/2);                                                                          % number of positive frequencies
                                    if mod(nbin,2)>0                                                                                % nbin is odd
                                        stftg(i,:,ilay)=[0 stftgl(i,2:nbinp) -stftgl(i,nbinp:-1:2) NaN(1,maxbinout-nbin)];          % force antisymmetry (odd nbin)
                                    else                                                                                            % nbin is even
                                        stftg(i,:,ilay)=[0 stftgl(i,2:nbinp-1) 0 -stftgl(i,nbinp-1:-1:2) NaN(1,maxbinout-nbin)];    % force antisymmetry (even nbin)
                                    end
                                end
                            end
                        else % real symmetric
                            if mod(lsym(ilay),2) % frequencies of this layer are shifted by 0.5 * frequency increment
                                for i=1:nfout                                                                                       % loop through frames to force conjugate symmetry
                                    nbin=metag(i,3);                                                                                % DFT length for this frame
                                    nbinp=1+floor(nbin/2);                                                                          % number of positive frequencies
                                    stftg(i,:,ilay)=[stftgl(i,1:nbin-nbinp+1) stftgl(i,nbinp-1:-1:1) NaN(1,maxbinout-nbin)];        % force symmetry
                                end
                            else % frequencies of this layer are aligned with STFT
                                for i=1:nfout                                                                                       % loop through frames to force conjugate symmetry
                                    nbin=metag(i,3);                                                                                % DFT length for this frame
                                    nbinp=1+floor(nbin/2);                                                                          % number of positive frequencies
                                    stftg(i,:,ilay)=[stftgl(i,1:nbinp) stftgl(i,nbin-nbinp+1:-1:2) NaN(1,maxbinout-nbin)];          % force symmetry
                                end
                            end
                        end
                    end

                end
            end                                                     % end of 1D or 2D interpolation: if strcmp(q.interpstft,'indep') ... else ... end
            % now sort out the interpolation of the metadata
            [xxi,xxf]=v_interval(taxout,taxin);                                         % i'th fixed frame centre, taxout(i), lies between taxin(xxi(i)) and taxin(xxi(i)+1)
            % interpolate the frame offset values
            switch q.interpof
                case 'lin'
                    metag(:,4)=meta(xxi,4).*(1-xxf)+meta(xxi+1,4).*xxf;                 % linearly interpolate the frame offset
            end
            % interpolate the frame scale factor values
            switch q.interpsc
                case 'lin'
                    metag(:,5)=meta(xxi,5).*(1-xxf)+meta(xxi+1,5).*xxf;                 % linearly interpolate the frame scale factor
                case 'log'
                    metag(:,5)=exp(log(meta(xxi,5)).*(1-xxf)+log(meta(xxi+1,5)).*xxf);  % linearly interpolate the log frame scale factor
            end
            % Calculate goup delay of output frames by using linear interpolation between the group delays of the input frames whose centres are either side of the
            % centre of the output frame while compensating for the starting sample of each of the frames.
            switch q.interpgd
                case 'lin'
                    metag(:,6)=(meta(xxi,1)+meta(xxi,6)).*(1-xxf)+(meta(xxi+1,1)+meta(xxi+1,6)).*xxf-metag(:,1);                 % linearly interpolate the group delay
                case 'linrep'
                    msk=xxf>0.5;                                                % mask for fixed frames closer to xxi(i)+1 than to xxi(i)
                    wtj=1-msk-(1-2*msk).*xxf;                                   % weight to apply to gdfj below: 1-xxf(i) if msk(i)=0 or xxf(i) if mask(i)=1
                    xxj=xxi+msk;                                                % taxout(i) is closest to  taxin(xxj(i))
                    xxk=xxi+1-msk;                                              % other end of interval
                    gdfj=v_modsym(meta(xxj,1)+meta(xxj,6),meta(xxj,3),taxout);  % add multple of DFT length to get assumed energy peak near centre of output frame
                    gdfk=v_modsym(meta(xxk,1)+meta(xxk,6),meta(xxk,3),gdfj);    % add multple of DFT length to get assumed energy peak near previous one
                    metag(:,6)=gdfj.*wtj+gdfk.*(1-wtj)-metag(:,1);              % group delay of fixed frame in samples: linear interpolate between gdfj and gdfk then compensate for start of frame
            end
            % compensate the stft values for the metadata and the starting sample of each frame
            for ilay=1:nlay
                if lsym(ilay)<0 % complex STFT so adjust group delay using metadata
                    for i=1:nfout % for now, apply the group delay frame by frame
                        nfft=metag(i,3);                     % DFT length of this frame
                        stftg(i,1,ilay)=stftg(i,1,ilay)-metag(i,4)*metag(i,3);   % subtract frame_offset*DFT_length from the DC value
                        stftg(i,1:nfft,ilay)=stftg(i,1:nfft,ilay).*exp(2i*pi/nfft*(interptzq*metag(i,1)+metag(i,6))*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1])/metag(i,5); % apply non-integer group delay (except to Nyquist frequency)
                    end
                end
            end
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