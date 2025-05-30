function [mfdr, cq, pa, naq, f0, h1h2, hrf] = extractVoiceFeatures(u, fs, gci)

%EXTRACTVOICEFEATURES extracts voice features from the glottal flow
%derivative  [mfdr, cq, pa, naq] = extractVoiceFeatures(u, fs, gci)
%
%  Inputs:  u        is the glottal flow in [cm^3/s]
%           fs       is the sampling frequency in Hz
%           gci      is the glottal closure instants indices (samples) of
%                    length Ng
%
% Outputs: 
% Features defined in Patel 2011 that are derived from glottal flow
% of length Ng
%
% MFDR  - Maximum Flow Declination Rate  [l/s^2]
% CQ    - Closed Quotient  [unitless]
% PA    - Pulse Amplitude  [cm^3/s]
% NAQ   - Normalized Amplitude Quotient [unitless]
%
% ___________________________
% not sure about these
% H1-H2  - Level difference between first and second harmonics
% Alpha  - Ratio between the summed energy between 50Hz-1kHz and 1-5 kHz
%          (computed from the long term average spectrum (LTAS))
%
% the rest of the parameters are obtained directly from sp
%
% leq   - Equivalent sound level
% Shimmer
% HNR  (harmonic-noise-ratio)
% Jitter
% MF0 (mean f0)


nfft = 2^12;

uu=[diff(u); 0]*fs;

mfdr = zeros(size(gci));
cq = zeros(size(gci));
pa = zeros(size(gci));
naq = zeros(size(gci));
opThres = 0.05;

tsrat=0.9;
terat=1-tsrat;
T=gci(2)-gci(1);

for ig=1:length(gci)
    nn=(-ceil(tsrat*T):ceil(terat*T))+gci(ig);  
    
    % If the indices exceed the signal boundaries then cut (and adjust T)
    nn=nn(nn>0); 
    nn=nn(nn<length(u));
    T=length(nn)-2;    % T should be 2 samples short of length nn
    
    
    uuseg=uu(nn);
    useg=u(nn);
    
    %useg = cumsum(uuseg)/fs;
    %useg=useg-linspace(0,useg(end),length(useg))';  % Enforce equal pressure at end
    %u(nn)=useg;
    
    % Maximum flow declination rate
    dpeak = -min(uuseg);  % flow derivative in cm^3/s^2
    % Maximum flow
    fac = max(useg) - min(useg);      % flow in cm^3/s
    % Pitch period 
    Ttime=T/fs;           % in seconds
    
    % Determining the duration of the open phase from the flow
    usegShift = useg-median(useg);
    
    chsign=diff(medfilt1(double(usegShift>opThres*max(usegShift)),7));  %Threshold the useg, median filter and find edges
    pch = find(chsign==1);    % Find where it goes to one (potential start of open phase)
    nch = find(chsign==-1);   % 
    for ii = 1:length(pch)  % This way we start looking at "open segments" where the first "opening" is
        nchIx = find(pch(ii)<nch);
        if ~any(nchIx)  % there is no closing at the end, so assume that it is the end of the segment
            segLen(ii) = length(useg) - pch(ii);
        else
            segLen(ii) = nch(nchIx(1)) - pch(ii);
        end;
    end;
    if isempty(pch)
        %warning('ba')
        cq(ig) = 0;   % Closed quotient is the 1 -  (the largest "open segment") / T
        mfdr(ig) = 0;   % in liters per second squared
        naq(ig) = 0;  % unit check: cm^3/s / (cm^3/s^2 * s) = unitless
        pa(ig) = 0;
        f0(ig) = 1/Ttime;
    else
        cq(ig) = 1 - max(segLen)/T;   % Closed quotient is the 1 -  (the largest "open segment") / T
        mfdr(ig) = dpeak/1000;   % in liters per second squared
        naq(ig) = fac/(dpeak*Ttime);  % unit check: cm^3/s / (cm^3/s^2 * s) = unitless
        pa(ig) = fac;
        f0(ig)=1/Ttime;
    end
    
    %Spectral parameters
    [pxx1,f] = pwelch(useg-mean(useg),hann(T),0,nfft,fs,'onesided');
    PSD1=10*log10(pxx1/max(pxx1));
    
    %     figure(1); clf;
    %     plot(PSD1);
    %     pause
    
    m=1;
    for k=f0(ig):f0(ig):3000
        BW=f0(ig)*.3;
        rangeF=(f>k-BW & f<k+BW);
        harm(m,1)= max(PSD1(rangeF));
        m=m+1;
    end
    h1h2(ig)=-harm(2,:);
    hrf(ig)=10*log10(sum(10.^(harm(2:end,:)/10)));

    if ig~=length(gci)
        T=gci(ig+1)-gci(ig);
    end
    %  Else, T stays the same (we assume that the last period is the same
    %  as the penultimate one.
    
    %T=gci(min(ig+1,length(gci)))-gci(ig);
end;
