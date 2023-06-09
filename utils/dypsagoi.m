function [gci,goi,gcic,goic,udash,crnmp] = dypsagoi(s,fs,opt)
%DYPSA   Derive glottal closure instances from speech [gci,goi] = (s,fs)
%   Note: Needs to be combined with a voiced-voiceless detector to eliminate
%   spurious closures in unvoiced and silent regions.
%
%   Inputs:
%   s     is the speech signal
%   fs    is the sampling frequncy
%   opt   (Optional) 'v' applies VUS detection
%
%   Outputs:
%   gci   is a vector of glottal closure sample numbers
%   gco   is a vector of glottal opening sample numbers derived from
%         an assumed constant closed-phase fraction
%
%   References:
%   [1]  P. A. Naylor, A. Kounoudes, J. Gudnason, and M. Brookes, �Estimation of Glottal Closure
%        Instants in Voiced Speech using the DYPSA Algorithm,� IEEE Trans on Speech and Audio
%        Processing, vol. 0, pp. 0�0, Jan. 2007 (accepted).
%   [2]  M. Brookes, P. A. Naylor, and J. Gudnason, �A Quantitative Assessment of Group Delay Methods
%        for Identifying Glottal Closures in Voiced Speech,� IEEE Trans on Speech & Audio Processing,
%        vol. 14, no. 2, pp. 456�466, Mar. 2006.
%   [3]  A. Kounoudes, P. A. Naylor, and M. Brookes, �The DYPSA algorithm for estimation of glottal
%        closure instants in voiced speech,� in Proc ICASSP 2002, vol. 1, Orlando, 2002, pp. 349�352.
%   [4]  C. Ma, Y. Kamp, and L. F. Willems, �A Frobenius norm approach to glottal closure detection
%        from the speech signal,� IEEE Trans. Speech Audio Processing, vol. 2, pp. 258�265, Apr. 1994.

% Algorithm Parameters
%       The following parameters are defined in voicebox()
%
%   dy_cpfrac=0.3;           % presumed closed phase fraction of larynx cycle
%   dy_cproj=0.2;            % cost of projected candidate
%   dy_cspurt=-0.45;         % cost of a talkspurt
%   dy_dopsp=1;              % Use phase slope projection (1) or not (0)?
%   dy_ewdly=0.0008;         % window delay for energy cost function term [~ energy peak delay from closure] (sec)
%   dy_ewlen=0.003;          % window length for energy cost function term (sec)
%   dy_ewtaper=0.001;        % taper length for energy cost function window (sec)
%   dy_fwlen=0.00045;        % window length used to smooth group delay (sec)
%   dy_fxmax=500;            % max larynx frequency (Hz) 
%   dy_fxmin=50;             % min larynx frequency (Hz) 
%   dy_fxminf=60;            % min larynx frequency (Hz) [used for Frobenius norm only]
%   dy_gwlen=0.0030;         % group delay evaluation window length (sec)
%   dy_lpcdur=0.020;         % lpc analysis frame length (sec)
%   dy_lpcn=2;               % lpc additional poles
%   dy_lpcnf=0.001;          % lpc poles per Hz (1/Hz)
%   dy_lpcstep=0.010;        % lpc analysis step (sec)
%   dy_nbest=5;              % Number of NBest paths to keep
%   dy_preemph=50;           % pre-emphasis filter frequency (Hz) (to avoid preemphasis, make this very large)
%   dy_spitch=0.2;           % scale factor for pitch deviation cost
%   dy_wener=0.3;            % DP energy weighting
%   dy_wpitch=0.5;           % DP pitch weighting
%   dy_wslope=0.1;           % DP group delay slope weighting
%   dy_wxcorr=0.8;           % DP cross correlation weighting
%   dy_xwlen=0.01;           % cross-correlation length for waveform similarity (sec)

%   Revision History: 
%
%   3.0 - 29 Jun 2006  - Rewrote DP function to improve speed
%   2.6 - 29 Jun 2006  - Tidied up algorithm parameters
%   2.4 - 10 Jun 2006  - Made into a single file aand put into VOICEBOX
%   2.3 - 18 Mar 2005  - Removed 4kHz filtering of phase-slope function 
%   2.2 - 05 Oct 2004  -  dpgci uses the slopes returned from xewgrdel
%                      -  gdwav from speech with fs<9000 is not filtered
%                      -  Various outputs and inputs of functions have been
%                         removed since now there is no plotting

%   Bugs:
%         1. Allow the projections only to extend to the end of the larynx cycle
%         2. Compensate for false pitch period cost at the start of a voicespurt
%         3. Should include energy and pahse-slope costs for the first closeure of a voicespurt
%         4. should delete candidates that are too close to the beginning or end of speech for the cost measures
%            currently this is 200 samples fixed in the main routine but it should adapt to window lengths of
%            cross-correlation, lpc and energy measures.
%         5. Should have an integrated voiced/voiceless detector
%         6. Allow dypsa to be called in chunks for a long speech file
%         7. Do forward & backward passes to allow for gradient descent and/or discriminative training
%         8. Glottal opening approximation does not really belong in this function
%         9. The cross correlation window is asymmetric (and overcomplex) for no very good reason
%        10. Incorporate -0.5 factor into dy_wxcorr and abolish erroneous (nx2-1)/(nx2-2) factor
%        11. Add in talkspurt cost at the beginning rather than the end of a spurt (more efficient)
%        12. Remove qmin>2 condition from voicespurt start detection (DYPSA 2 compatibility) in two places
%        13. Include energy and phase-slope costs at the start of a voicespurt
%        14. Single-closure voicespurt are only allowed if nbest=1 (should always be forbidden)
%        15. Penultimate closure candidate is always acceptd
%        16. Final element of gcic, Cfn and Ch is unused
%        17. Needs to cope better with irregular voicing (e.g. creaky voice)
%        18. Should give non-integer GCI positions for improved accuracy
%        19. Remove constraint that first voicespurt cannot begin until qrmax after the first candidate

%      Copyright (C) Tasos Kounoudes, Jon Gudnason, Patrick Naylor and Mike Brookes 2006
%      Version: $Id: dypsa.m,v 3.3 2006/07/28 07:40:57 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<3)
    opt='';
end
if(nargin>=3)
    if ~(opt=='v')
        error('Argument ''opt'' can only be ''v'' for now. Stay tuned for more exciting options.');
    end
end

% Extract algorithm constants from VOICEBOX
dy_preemph=voicebox('dy_preemph');
dy_lpcstep=voicebox('dy_lpcstep');
dy_lpcdur=voicebox('dy_lpcdur');
dy_dopsp=voicebox('dy_dopsp');              % Use phase slope projection (1) or not (0)?
dy_ewtaper=voicebox('dy_ewtaper');        % Prediction order of FrobNorm method  in seconds
dy_ewlen=voicebox('dy_ewlen');        % windowlength of FrobNorm method  in seconds
dy_ewdly=voicebox('dy_ewdly');        % shift for assymetric speech shape at start of voiced cycle
dy_cpfrac=voicebox('dy_cpfrac');        % presumed ratio of larynx cycle that is closed
dy_lpcnf=voicebox('dy_lpcnf');          % lpc poles per Hz (1/Hz)
dy_lpcn=voicebox('dy_lpcn');            % lpc additional poles

lpcord=ceil(fs*dy_lpcnf+dy_lpcn);       % lpc poles

spp = s;

% DC removal by high-pass filter - REINTRODUCED %
[b a] = butter(2,10/(fs/2),'high');
spf = filter(b,a,spp);

% MRPT's parameters
IF = 'iaif';
IF = 'mark';
postGOI = 1;

if(strcmp(IF,'mark'))
    %load('/Users/mark/Work/Project/Excitation Modelling/ddm_aplawd/models/protMelcepstNormudash1.mat');
    %protudash = 0.5*(protudash(1:400)+protudash(401:end));
    %iprotudash = lsinvfilt(protudash,400,0);
    
    iprotudash=[
        -0.533002849324892; 
        0.622127983250758;
        -0.058083131786265;
        -0.025840951313764;
        -0.002462576708046;
        0.002487699174393;
        -0.000280997739956;
        -0.002140389900568;
        -0.001615864970452;
        -0.000987664078621;
        -0.001220503224164;
        -0.001438863736837;
        -0.001176551923861;
        -0.000537854598417;
        0.000010322336443;
        0.000105351146253];
    
    s_used = fftfilt(iprotudash,spf);
    
    % perform LPC analysis, AC method with Hamming windowing
    [ar, e, Ts] = lpcauto(s_used,lpcord,floor([dy_lpcstep dy_lpcdur]*fs));
    udash = lpcifilt(spf,ar,Ts);    % Pad gives same alignment as original LPC residual.
elseif(strcmp(IF,'iaif'))
    s_used=filter([1 -exp(-2*pi*dy_preemph/fs)],1,s);
    udash = myiaif(s,fs,20,4,20,1);
end

% Calculate SWT
nlev=3;
nu = length(udash);
nU=(2^nlev)*ceil(nu./(2^nlev));
[Lo_D Hi_D] = wfilters('bior1.5','d');    
[swa swd] = swtalign([udash; zeros(nU-nu,1)],nlev,Lo_D,Hi_D);
swa=swa(:,1:nu); swd=swd(:,1:nu);

% Find multiscale products
mp = prod(swd)';

% Find third roots
nmp = mp;
nmp(find(nmp>0)) = 0;   % Half-wave rectify on negative half of mp for GCI
crnmp = nthroot(nmp,3);
r=crnmp;

% compute the group delay function:  EW method from reference [2] above
[zcr_cand,sew,gdwav,toff]=xewgrdel(r,fs); 
gdwav=-[zeros(toff,1); gdwav(1:end-toff)];
zcr_cand=[round(zcr_cand), ones(size(zcr_cand))];   %flag zero crossing candidates with ones

sew=0.5+sew';  %the phase slope cost of each candidate

pro_cand=[];
if dy_dopsp ~= 0
    pro_cand = psp(gdwav,fs);
    pro_cand = [pro_cand, zeros(length(pro_cand),1)]; %flag projected candidates with zeros
    sew =      [sew zeros(1,size(pro_cand,1))];      %the phase slope cost of a projected candidate is zero
end;

%Sort the zero crossing and projected candidates together and remove any candidates that
%are within 200 samples [0.01 s] from the speech boundary
[gcic,sin] = sortrows([zcr_cand; pro_cand],1);  
sew=sew(sin);
sin=find(and(0.01*fs<gcic,gcic<length(gdwav)-0.01*fs));
gcic=gcic(sin,:);
sew=sew(sin);

% compute the frobenious norm function used for a cost in the DP
fnwav=frobfun(s_used,dy_ewtaper*fs,dy_ewlen*fs,dy_ewdly*fs);    % Original

% Closed phase detection
u = filter(1,[1, -0.99],udash);     % Don't know why this works - better behaved than udash?
aencost = zeros(length(gcic),1);
for i=1:length(gcic)-1
    tmp = u(gcic(i,1):gcic(i+1,1));
    aencost(i)=(mean(tmp));
end
aencost = 0.5*aencost/mean(abs(aencost));
cencost = [0; aencost(1:end-1)];

[gci] = dpgci(gcic, udash(:), sew, fnwav, fs, aencost, double(opt=='v'));

% Remove GCIs from candidate set
if(nargout>1)
    [goic I] = setdiff(gcic(:,1),gci(:));
    goic = [goic gcic(I,2)];
    sewo = sew(I);
    cencost = cencost(I);
end

% Refine GCIs as per hybrid method.
k=v_findpeaks(-crnmp,'',50);
k=k(:);
gci=gci(:);
tol=10; % Was 15 - not sure how to optimise this.
[i j] = find(abs(repmat(k,1,length(gci)) - repmat(gci',length(k),1))<tol);
j2 = setdiff(1:length(gci),j);
gci=[k(i); gci(j2)];
gci = sort(gci)';

% MRPT: GOI Postprocessing
if(nargout>1)
    %[goi] = dpgci(goic, udash(:), sewo, fnwav,fs);
    [goi] = dpgci(goic, udash(:), sewo, fnwav,fs,cencost);
    
    % NEEDS FIXING - DOESN'T ALWAYS GIVE EQUAL GCIS AS GOIS
    if(postGOI)
        %for i=1:2   % New and temporary
        % When a GOI doesn't follow a GCI, add at previous period and remove any
        % stray GOIs (should run GCI-GOI-GCI-GOI...).
        % No VUS required: if run on GCI then following code sorts it.
        k = [ones(size(gci)) -ones(size(goi))];
        x = [[gci goi];k];

        [gi I] = sort(x(1,:),2);
        k = k(I);
        x = x(:,I);
        p = find(fftfilt([1 1],k)>0)-1; % Those to add
        p(find(p==0))=1;
        n = find(fftfilt([1 1],k)<0);   % Those to remove

        x(2,n) = 0;     % Flag all extras for removal
        x = [x zeros(2,length(p))];
        for i=1:length(p)
            [a I] = find(x(2,1:p(i))==-1,1,'last'); % Find closest previous GOI
            if(I>1)
                x(:,end+i) = [x(1,p(i))+(x(1,I)-(x(1,I-1)));-1];
            else
                x(:,end+i) = [x(1,p(i))+(x(1,I+1)-x(1,I));-1];
            end
        end
        [f_ I] = find(x(2,:)==-1);
        goi = sort(x(1,I));
        %end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = psp(g,fs)
%PSP  Calculates the phase slope projections of the group delay function
%   Z = PSP(G) computes the 

%   Author(s): P. A. Naylor
%   Copyright 2002 Imperial College of Science Technology and Medicine, London
%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

g = g(:);

gdot = [diff(g);0];
gdotdot = [diff(gdot);0];

% find the turning points  as follows: [tp_number, index_of_tp, min(1) or max(-1), g(index_of_tp)]
turningPoints = zcr(gdot);

turningPoints = [[1:length(turningPoints)]', turningPoints, sign(gdotdot(turningPoints)), g(turningPoints)];

% useful for debug/plotting
%tplot = zeros(length(g),1);
%tplot(turningPoints(:,1)) = turningPoints(:,2);

% find any maxima which are < 0 
negmaxima = turningPoints(find(turningPoints(:,3) == -1 & turningPoints(:,4) < 0 & turningPoints(:,1)~=1),:);  %Change 01.05.2003 JG: The first row can't be included

% find the midpoint between the preceding min and the negative max
nmi = negmaxima(:,1);
midPointIndex = turningPoints(nmi-1,2) + round(0.5*(turningPoints(nmi,2) - turningPoints(nmi-1,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
nz = midPointIndex - round(midPointValue);

% find any minima which are > 0 
posminima = turningPoints(find(turningPoints(:,3) == 1 & turningPoints(:,4) > 0),:);

% find the midpoint between the positive min and the following max
pmi = posminima(:,1); 

%Remove last midpoint if it is the last sample
if ~isempty(pmi), if pmi(end)==size(turningPoints,1), pmi=pmi(1:end-1); end; end;

midPointIndex = turningPoints(pmi,2) + round(0.5*(turningPoints(pmi+1,2) - turningPoints(pmi,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
pz = midPointIndex - round(midPointValue);

z = sort([nz;pz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function i = zcr(x, p)
%ZCR  Finds the indices in a vector to  zero crossings
%   I = ZCR(X) finds the indices of vector X which are closest to zero-crossings.
%   I = ZCR(X, P) finds indices for positive-going zeros-crossings for P=1 and
%   negative-going zero-crossings for P=0.

x = x(:);

if (nargin==2)
    if (p==0) 
        z1 = zcrp(x);   % find positive going zero-crossings
    elseif (p==1) 
        z1 = zcrp(-x);  % find negative going zero-crossings
    else
        error('ZCR: invalid input parameter 2: must be 0 or 1');
    end
else
    z1 = [zcrp(x); zcrp(-x)];
end

% find crossings when x==0 exactly
z0 = find( (x(1:length(x))==0) & ([x(2:length(x));0] ~= 0));

% concatenate and sort the two types of zero-crossings
i = sort([z0; z1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zz = zcrp(xx)  %only used in zcr
% find positive-going zero-crossing
z1 = find(diff(sign(xx)) == -2);
% find which out of current sample or next sample is closer to zero
[m, z2] = min([abs(xx(z1)), abs(xx(z1+1))], [], 2);
zz =  z1 -1 + z2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frob]=frobfun(sp,p,m,offset)

% [frob]=frobfun(sp,p,m)
% 
% sp is the speech signal assumed to be preemphasised
% p  is the prediction order  : recomended to be 1 ms in above paper
% m  is the window length     : recomended to be 1 ms in above paper
% offset is shift for assymetric speech shape at start of voiced cycle -
% default 1.5ms.
%
% This function implements the frobenius norm based measure C defined in [4] below.
% It equals the square of the Frobenius norm of the m by p+1 data matrix divided by p+1
%
% Reference:
%   [4]  C. Ma, Y. Kamp, and L. F. Willems, �A Frobenius norm approach to glottal closure detection
%        from the speech signal,� IEEE Trans. Speech Audio Processing, vol. 2, pp. 258�265, Apr. 1994.


%   Author(s): J. Gudnason, P
%   Copyright 2002 Imperial College of Science Technology and Medicine, London
%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

%force p m and offset to be integers
p=round(p);
m=round(m);
offset=round(offset);

w=(p+1)*ones(1,m+p);
w(1:p)=1:p;
w(m+1:p+m)=p:-1:1;

w=w./(p+1); 
frob=filter(w,1,sp.^2);
frob(1:(round((p+m-1)/2) + offset))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function goi=simplegci2goi(gci,pr)

% Estimate glottal opening instants by assuming a fixed closed-phase fraction

gci=round(gci);
maxpitch=max(medfilt1(diff(gci),7));

% calculate opening instants
for kg=1:length(gci)-1
    goi(kg)=gci(kg)+min(pr*(gci(kg+1)-gci(kg)),pr*maxpitch);
end;
kg=kg+1;
goi(kg)=round(gci(kg)+pr*(gci(kg)-gci(kg-1)));  %use the previous pitch period instead
goi=round(goi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tew,sew,y,toff]=xewgrdel(u,fs)

% implement EW group delay epoch extraction

%dy_gwlen=voicebox('dy_gwlen');          % group delay evaluation window length
%dy_fwlen=voicebox('dy_fwlen');          % window length used to smooth group delay

dy_gwlen = 0.002;      % 0.002 works well
dy_fwlen = 0.00045;     % 0.00045 works well

% perform group delay calculation

gw=2*floor(dy_gwlen*fs/2)+1;            % force window length to be odd
ghw=window('hamming',gw,'s');
ghw = ghw(:);                           % force to be a column (dmb thinks window gives a row - and he should know as he wrote it!)
ghwn=ghw'.*(gw-1:-2:1-gw)/2;            % weighted window: zero in middle

u2=u.^2;
yn=filter(ghwn,1,u2);
yd=filter(ghw,1,u2);
yd(abs(yd)<eps)=10*eps;                 % prevent infinities
y=yn(gw:end)./yd(gw:end);               % delete filter startup transient
toff=(gw-1)/2;
fw=2*floor(dy_fwlen*fs/2)+1;            % force window length to be odd
if fw>1
    daw=window('hamming',fw,'s');
    y=filter(daw,1,y)/sum(daw);         % low pass filter 
    toff=toff-(fw-1)/2;
end
[tew,sew]=zerocros(y,'n');              % find zero crossings

tew=tew+toff;                           % compensate for filter delay and frame advance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cfn=fnrg(gcic,frob,fs)

%Frobenious Energy Cost

dy_fxminf=voicebox('dy_fxminf');
frob=frob(:)';
mm=round(fs/dy_fxminf);
mfrob=maxfilt(frob,1,mm);
mfrob=[mfrob(floor(mm/2)+1:end) max(frob(end-ceil(mm/2):end))*ones(1,floor(mm/2))];
rfr=frob./mfrob;
Cfn=0.5-rfr(round(gcic));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gci=dpgci(gcic, s, Ch, fnwav, fs, encost, vus)

%DPGCI   Choose the best Glottal Closure Instances with Dynamic Programming
%   gci=dpgci(gcic, s(:), Ch, fnwav, fs) returns vectors of sample indices corresponding
%   to the instants of glottal closure in the speech signal s at sampling frequency fs Hz.
%
%   Inputs:
%   gcic    is a matrix whos first column are the glottal closure instance candidates and
%           the second column is 1 if the corresponding gci is derived from a zero crossing 
%           but zero if the gci is from a a projected zero crossing
%   s       is the speech signal - MUST be a column vector
%   Ch      the phase slope cost of every candidate
%   fnwav   is the frobenious norm function of s
%   fs      is the sampling frequncy
%   encost  MRPT: cost for anticausal energy (for GCIs) and causal energy
%           (GOIs)
%
%   Outputs:
%   gci     is a vector of glottal closure instances chosen by the DP



%   Author(s): J. Gudnason, P.A. Naylor and D.M. Brookes
%   Copyright 2003 Imperial College London
%   Revision History: 
%   Bugs:  Constants are hardwired but defined in a structure like pv (defined in grpdelpv)
%         

% get algorithm parameters from voicebox()

dy_fxmin=voicebox('dy_fxmin');        % min larynx frequency (Hz)
dy_fxmax=voicebox('dy_fxmax');        % min larynx frequency (Hz)
dy_xwlen=voicebox('dy_xwlen');        % cross-correlation length for waveform similarity (sec)
dy_nbest=voicebox('dy_nbest');        % Number of NBest paths to keep
dy_spitch=voicebox('dy_spitch');              % scale factor for pitch deviation cost
wproj=voicebox('dy_cproj');           % cost of projected candidate
dy_cspurt=voicebox('dy_cspurt');           % cost of a talkspurt
dy_wpitch=voicebox('dy_wpitch');           % DP pitch weighting
dy_wener=voicebox('dy_wener');           % DP energy weighting
dy_wslope=voicebox('dy_wslope');           % DP group delay slope weighting
dy_wxcorr=voicebox('dy_wxcorr');           % DP cross correlation weighting

%Constants
Ncand=length(gcic);
sv2i=-(2*dy_spitch^2)^(-1);              % scale factor for pitch deviation cost

%Limit the search:
qrmin=ceil(fs/dy_fxmax);
qrmax=floor(fs/dy_fxmin);

% % MRPT: Cost saving
% savecosts = 1;
if(nargin<7)
    vus =0;
end

%Cost and tracking r = current, q = previous, p = preprevious
cost=zeros(Ncand, dy_nbest); cost(:,:)=inf;    %Cost matrix, one row for each candidate
maxcost=zeros(Ncand,1); maxcost(:,:)=inf;   %Maximum cost in each row
imaxcost=ones(Ncand,1);                     %Index of maximum cost

prev = ones(Ncand, dy_nbest);                  %index of previous, q candidates
ind = ones(Ncand, dy_nbest);                   %index of p in row q (from prev)
qbest = [zeros(Ncand,1), ones(Ncand,2)]; % the minimum cost in any previous q [cost,q,i]

% Temp: replace fnwav with crnmp
Cfn=fnrg(gcic(:,1),fnwav,fs);  %Frob.Energy Cost

%Add start and end state
% === should probably delete candidates that are too close to either end of the input
% === why do we ever need the additional one at the tail end ?
gcic=[[gcic(1,1)-qrmax-2 0];gcic;[gcic(end,1)+qrmax+2 0]];
Cfn=[0 Cfn 0];
Ch = [0 Ch 0];

% first do parallelized version

nxc=ceil(dy_xwlen*fs);       % cross correlation window length in samples
% === should delete any gci's that are within this of the end.
% === and for energy window
% rather complicated window specification is for compatibility with DYPSA 2
% === +1 below is for compatibility - probably a bug
wavix=(-floor(nxc/2):floor(nxc/2)+1)';                 % indexes for segments [nx2,1]
nx2=length(wavix);
sqnx2=sqrt(nx2);

if(nargin<6)
    g_cr=dy_wener*Cfn+dy_wslope*Ch+wproj*(1-gcic(:,2))';  % fixed costs
else
    dy_wamp = 0.5;
    g_cr=dy_wener*Cfn+dy_wslope*Ch+wproj*(1-gcic(:,2))'+dy_wamp*[0 encost' 0];  % fixed costs
end


g_n=gcic(:,1)';                  % gci sample number [1,Ncand+2]
g_pr=gcic(:,2)';                 % unprojected flag [1,Ncand+2]
g_sqm=zeros(1,Ncand+1);         % stores: sqrt(nx2) * mean for speech similarity waveform
g_sd=zeros(1,Ncand+1);         % stores: 1/(Std deviation * sqrt(nx2)) for speech similarity waveform
f_pq=zeros((Ncand+1)*dy_nbest,1);   % (q-p) period for each node
f_c=repmat(Inf,(Ncand+1)*dy_nbest,1);    % cumulative cost for each node - initialise to inf
f_c(1)=0;                       % initial cost of zero for starting node
% f_costs=zeros(Ncand*dy_nbest,6);   % === debugging only remember costs of candidate
f_f=ones((Ncand+1)*dy_nbest,1);    % previous node in path
f_fb=ones((Ncand+1),1);    % points back to best end-of-spurt node
fbestc=0;                       % cost of best end-of-spurt node

qmin=2;
for r=2:Ncand+1   
%     if r==86
%         r;
%     end
    r_n=g_n(r);             % sample number of r = current candidate
    rix=dy_nbest*(r-1)+(1:dy_nbest);    % index range within node variables
    
    % determine the range of feasible q candidates
    qmin0=qmin;
    qmin=find(g_n(qmin0-1:r-1)<r_n-qrmax);      % qmin is the nearest candidate that is >qrmax away
    qmin=qmin(end)+qmin0-1;             % convert to absolute index of first viable candidate
    qmax=find(g_n(qmin-1:r-1)<=r_n-qrmin);      % qmax is the nearest candidate that is >=qrmin away
    qmax=qmax(end)+qmin-2;
    
    
    % calculate waveform similarity cost measure statistics
    sr=s(r_n+wavix);        % note s MUST be a column vector so sr is also
    wsum=sum(sr);
    g_sqm(r)=wsum/sqnx2;                % mean * sqrt(nx2)
    g_sd(r)=1/sqrt(sr.'*sr-wsum^2/nx2);   % 1/(Std deviation * sqrt(nx2))
    
    % now process the candidates
    
    if qmin<=qmax
        qix=qmin:qmax;      % q index
        nq=length(qix);
        % === should integrate the -0.5 into dy_wxcorr
        % === the factor (nx2-1)/(nx2-2) is to compensate for a bug in swsc()
        q_cas=-0.5*(nx2-1)/(nx2-2)*dy_wxcorr*(sum(s(repmat(g_n(qix),nx2,1)+repmat(wavix,1,nq)).*repmat(sr,1,nq),1)-g_sqm(qix)*g_sqm(r)).*g_sd(qix)*g_sd(r);
        % compare: i=35; Ca=swsc(g_n(qix(i)),g_n(r),s,fs); [i qix(i) r  g_n(qix(i)) g_n(r) dy_wxcorr*Ca q_cas(i)]
        
        % now calculate pitch deviation cost
        
        fix = 1+(qmin-1)*dy_nbest:qmax*dy_nbest;    % node index range
        f_qr=repmat(r_n-g_n(qix),dy_nbest,1);    % (r-p) period for each node
        f_pr=f_qr(:)+f_pq(fix);
        % === could absorb the 2 into sv2i
        f_nx=2-2*f_pr./(f_pr+abs(f_qr(:)-f_pq(fix)));
        f_cp=dy_wpitch*(0.5-exp(sv2i*f_nx.^2));
        % === fudge to match dypsa2.4 - could more efficiently be added
        % === onto the cost of a talkspurt end
        % === should be a voicebox parameter anyway
        f_cp(f_pq(fix)==0)=dy_cspurt*dy_wpitch;
        
        % now find the N-best paths
        
        [r_cnb,nbix]=sort(f_c(fix)+f_cp+reshape(repmat(q_cas,dy_nbest,1),nq*dy_nbest,1));
        f_c(rix)=r_cnb(1:dy_nbest)+g_cr(r);     % costs
        f_f(rix)=nbix(1:dy_nbest)+(qmin-1)*dy_nbest;       % traceback nodes
        f_pq(rix)=f_qr(nbix(1:dy_nbest));       % previous period
        % === f_costs is only for debugging
%         r;
%         f_costs(rix,1)=f_c(fix(nbix(1:dy_nbest)));
%         f_costs(rix,2)=wproj*(1-gcic(r,2));
%         f_costs(rix,3)=f_cp(nbix(1:dy_nbest));
%         f_costs(rix,4)=dy_wener*Cfn(r);
%         f_costs(rix,5)=dy_wslope*Ch(r);
        if(vus)   % Currently no weights
%             f_costs(rix,1)=reshape(q_cas(1+floor((nbix(1:dy_nbest)-1)/dy_nbest)),dy_nbest,1)./dy_wxcorr;    % Waveform similarity
%             f_costs(rix,2)=f_cp(nbix(1:dy_nbest))./dy_wpitch;                                               % Pitch deviation
%             f_costs(rix,3)=(1-gcic(r,2));      % Project candidate cost - weight in paper is 0.4 but weight is 0.2 as (1-gcic(r,2)) is bounded in [0,1].
%             f_costs(rix,4)=Cfn(r);                                                                          % Frobenius Energy cost
%             f_costs(rix,5)=Ch(r);                                                                           % Phase slope deviation cost

            f_costs(rix,1)=reshape(q_cas(1+floor((nbix(1:dy_nbest)-1)/dy_nbest))./dy_wxcorr,dy_nbest,1);    % Waveform similarity
            f_costs(rix,2)=f_cp(nbix(1:dy_nbest))./dy_wpitch;      % Pitch deviation
            f_costs(rix,3)=(1-gcic(r,2))/2;         % Project candidate cost - weight in paper is 0.4 but weight is 0.2 as (1-gcic(r,2)) is bounded in [0,1].
            f_costs(rix,4)=Cfn(r);             % Frobenius Energy cost
            f_costs(rix,5)=Ch(r);             % Phase slope deviation cost
        end

        % check cost of using this candidate as the start of a new spurt
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again below)
        iNb=rix(end);        
        if (qmin>2) && (f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2))<f_c(iNb))        % compare with worst of Nbest paths
            f_f(iNb)=f_fb(qmin-1);
            % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
            % === this is probably a bug
            f_c(iNb)=f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2));     % replace worst of the costs with start voicespurt cost
            f_pq(iNb)=0;                    % false pq period
        end
        if f_c(rix(1))<fbestc
            f_fb(r)=rix(1);                          % points to the node with lowest end-of-spurt cost
            % === should compensate for the pitch period cost incurred at the start of the next spurt
            % === note that a node can never be a one-node voicespurt on its own unless dy_nbest=1
            % since the start voices[purt option replaced the worst Nbest cost. This is probably good but
            % is a bit inconsistent.
            fbestc=f_c(rix(1));
        else
            f_fb(r)=f_fb(r-1);
        end
    else            % no viable candidates - must be the start of a voicespurt if anything
        % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
        % === this is probably a bug
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again above)
        if (qmin>2)
            f_c(rix(1))=f_c(f_fb(qmin-1))+wproj*(1-gcic(r,2));  % cost of new voicespurt
            f_f(rix)=f_fb(qmin-1);                              % traceback to previous talkspurt end
            f_pq(rix)=0;                                        % previous period
        end
        f_fb(r)=f_fb(r-1);                                  % cannot be the end of a voicespurt
    end
end

% now do the traceback

gci = zeros(1,Ncand+1);

if(vus)
    mycost = zeros(Ncand+1,5);
end

% === for compatibility with dypsa2, we force the penultimate candidate to be accepted
% === should be: i=f_fb(Ncand+1) but instead we pick the best of the penultimate candidate
i=rix(1)-dy_nbest;
if f_c(i-dy_nbest+1)<f_c(i)     % check if start of a talkspurt
    i=i-dy_nbest+1;
end
k=1;
while i>1
    j=1+floor((i-1)/dy_nbest);          % convert node number to candidate number
    gci(k)=g_n(j);
    if(vus)
        mycost(k,:) = f_costs(i,:);
    end
    i=f_f(i);
    k=k+1;
end
gci=gci(k-1:-1:1);           % put into ascending order

% VUS detector
if(vus)
    mycost=mycost(k-1:-1:1,:);

    % Smooth waveform similarity
    dy_fwlen = 0.001;
    fw=2*floor(dy_fwlen*fs/2)+1;            % force window length to be odd
    daw=window('hamming',fw,'s');
    y=filter(daw,1,mycost(:,1))/sum(daw);         % low pass filter
    delay = (length(daw)+1)/2;
    y = [y(delay+1:end); zeros(delay, 1)];

    iV = y<-0.3;
%     figure;
%     myax(1)=subplot(2,1,1);plot(s/max(abs(s)));hold on;plot(gci,mycost);legend('Speech', 'Waveform sim.', 'Pitch dev.', 'Proj. cand', 'Frob. ener.', 'PS dev.');stem(gci,ones(size(gci)),'g');title('u''(n), DYPSA GCIs and costs');
%     myax(2)=subplot(2,1,2);stem(gci,iV);%hold on;stem(evalin('base','sgci'),-ones(size(evalin('base','sgci'))),'r');title('GMM-chosen GCIs and reference');
%     linkaxes(myax,'x');
    
    % Kill any single or double bursts. Leaves holes in bursts length 3.
    iVr = [zeros(2,1); iV(1:end-2)] | [iV(2+1:end); zeros(2,1)];
    % Insert missing singles in both iV and (iV&iVr)
    iVi = (fftfilt([1 0 1], double((iV&iVr | iV)==2)));
    % Filter out
    %gci = gci(([iVi(2:end); 0] | (iV&iVr))==1);
    gciv = gci(([iVi(2:end); 0] | (iV&iVr)));
    gciu = setdiff(gci,gciv);
    
    gci = gciv;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = swtalign(x,n,varargin)
%SWT Discrete stationary wavelet transform 1-D. Differs from swt(...) in
%alignment of multiscale coefficients

nbIn = nargin;
if nbIn < 3
  error('Not enough input arguments.');
elseif nbIn > 4
  error('Too many input arguments.');
end
if errargt(mfilename,n,'int'), error('*'), end

% Use row vector.
x = x(:)';
s = length(x);
pow = 2^n;
if rem(s,pow)>0
    sOK = ceil(s/pow)*pow;
    msg = strvcat(...
            ['The level of decomposition ' int2str(n)],...
            ['and the length of the signal ' int2str(s)],...
            'are not compatible.',...
            ['Suggested length: ' int2str(sOK)],...
            '(see Signal Extension Tool)', ...
            ' ', ...
            ['2^Level has to divide the length of the signal.'] ...
            );
    errargt(mfilename,msg,'msg');
    varargout = {[] };
    return
end

% Compute decomposition filters.
if nargin==3
    [lo,hi] = wfilters(varargin{1},'d');
else
    lo = varargin{1};   hi = varargin{2};
end

% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

% Compute stationary wavelet coefficients.
evenoddVal = 0;
evenLEN    = 1;
swd = zeros(n,s);
swa = zeros(n,s);
for k = 1:n

    % Extension.
    lf = length(lo);
    x  = wextend('1D',modeDWT,x,lf/2);

    % Decomposition.
    swd(k,:) = wkeep1(wconv1(x,hi),s,lf);   % Default last arg was lf+1.
    swa(k,:) = wkeep1(wconv1(x,lo),s,lf);
    
    % upsample filters.
    lo = dyadup(lo,evenoddVal,evenLEN);
    hi = dyadup(hi,evenoddVal,evenLEN);

    % New value of x.
    x = swa(k,:);

end

if nargout==1
    varargout{1} = [swd ; swa(n,:)];
elseif nargout==2
    varargout = {swa,swd};
end     

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [udash ar Ts u] = myiaif(sp,fs,p,g,r,h)

%   Implementation of the Iterative Adaptive Inverse Filtering (IAIF)
%   Algorithm for Glottal Wave Analysis
%
%   [udash ar Ts u] = iaif(sp,fs,p,g,r,h)
%
%   Inputs:
%       sp      Nx1 vector speech signal
%       fs      Sampling freq (Hz)
%       p       (Optional) First vocal tract LPC order (8-12, default=10)
%       g       (Optional) Glottal source LPC order (2-4, default=4)
%       r       (Optional) Second vocal tract LPC order (8-12, default=10)
%       h       (Optional) 0-don't HPF at 60 kHz (default=1)
%
%   Outputs:
%       udash   Nx1 vector of glottal flow derivative
%       ar      
%       Ts
%       u       Nx1 vector of glottal flow
%
%   Notes:
%       1. Based upon algorithm definition in Alku1999.
%       2. LPC orders specified for 8k sampling only; inputs resampled if
%          necessary.
%       3. 'Integrator' assumed to mean 'slightly leaky integrator'.
%       4. LPC overlap of 50% assumed.
%       5. Highpass order assumed to be 1024.
%       6. Mark's preferred parameter sets: 8k:   p=8;g=2;r=8;
%                                           20k:  p=20;g=4;r=20;
%
%   External Functions:
%       Functions lpcauto and lpcifilt in Voicebox required.
%
%   References:
%       P. Alku, "Glottal Wave Analysis with Pitch Synchronous Iterative
%       Adaptive Filtering," Speech Communication, 1992, 11(2-3), 109-118.
%       
%       P. Alku, H. Tiitinen and R. Naatanen, "A Method for Generating
%       Natural-Sounding Speech Stimuli for Cognitive Brain Research,"
%       Clinical Neurophysiology, May 1999, 110(8), 1329-1333.
%
%**************************************************************************
% Author:           M. R. P. Thomas 
% Date:             28 April 2009
% Last Modified:    29 April 2009
%**************************************************************************

lpcdur=0.032;       % As in Alku1992a. Not stated in Alku1999.
lpcstep=lpcdur/2;   % 50% overlap assumed.

%if~(fs==8000)
%    warning('fs!=8000. Consider resampling');
%end

% LPC Orders. 8k recommendation: 8<=p<=12, 2<=g<=4, 8<=r<=12.
% Alku1992a set
if(nargin<6)
    h = 1;
end
if(nargin<5)
    r=10;
end
if(nargin<4)
    g=4;
end
if(nargin<3)
    p=10;
end

% MRPT's Alternative LPC parameter sets. Can it be made fs-dependent?
% 8k:   p=8;g=2;r=8;
% 20k:  p=25;g=4;r=25;

% 'Integrator' undefined - assume slightly leaky.
intb=1;
inta=[1 -0.95];

% 1. Highpass - 60k in Alku1999 and 30k in PSIAIF. May even make external.
if(h==1)
    b = fir1(1024,60/(fs/2),'high');
    spf=fftfilt(b,sp);
    delay=round(mean(grpdelay(b)));
    spf = [spf(delay+1:end); zeros(delay, 1)];
else
    spf=sp;
end

% 2. Order-1 LPC
[Hg1, eHg1, TsHg1] = lpcauto(spf,1,floor([lpcstep lpcdur]*fs));

% 3. Inverse-filtering
spHg1 = lpcifilt(spf,Hg1,TsHg1);

% 4. p-th-order LPC
[Hvt1, eHvt1, TsHvt1] = lpcauto(spHg1,p,floor([lpcstep lpcdur]*fs));

% 5. Inverse-filtering
spHvt1 = lpcifilt(spf,Hvt1,TsHvt1);

% 6. Integrate
spHvt1_I = filter(intb,inta,spHvt1);

% 7. g-th order LPC
[Hg2, eHg2, TsHg2] = lpcauto(spHvt1_I,g,floor([lpcstep lpcdur]*fs));

% 8. Inverse-filtering
spHg2 = lpcifilt(spf,Hg2,TsHg2);

% 9. Integrate - not done in Alku1992 but is in Alku1999
spHg2_I = filter(intb,inta,spHg2);

% 10. r-th order LPC
[Hvt2, eHvt2, TsHvt2] = lpcauto(spHg2_I,r,floor([lpcstep lpcdur]*fs));

% 11. Inverse-filtering
spHvt2 = lpcifilt(spf,Hvt2,TsHvt2);
udash = spHvt2;

ar=Hvt2;
Ts = TsHvt2;

% 12. Integration
g = filter(intb,inta,spHvt2);
u = g;


