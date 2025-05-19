function [c,tc,y,mc]=w_melcepstx(s,fs,w,nc,p,n,inc,fl,fh)
% enhanced version of v_melcepst with y and mc outputs and 'S' and 'f' options
%V_MELCEPST Calculate the mel cepstrum of a signal C=(S,FS,W,NC,P,N,INC,FL,FH)
%
%
% Simple use: (1) c=v_melcepst(s,fs)          % calculate mel cepstrum with 12 coefs, 256 sample frames
%			  (2) c=v_melcepst(s,fs,'E0dD')   % include log energy, 0th cepstral coef, delta and delta-delta coefs
%
% Inputs:
%     s	  speech signal (or stft(time,freq) with 'S' option)
%     fs  sample rate in Hz (default 11025)
%     w   mode string (see below)
%     nc  number of cepstral coefficients excluding 0'th coefficient [default 12]
%     p   number of filters in filterbank [default: floor(3*log(fs)) =  approx 2.1 per ocatave]
%     n   length of frame in samples [default power of 2 < (0.03*fs)]
%     inc frame increment in samples [default n/2]
%     fl  low end of the lowest filter as a fraction of fs [default = 0]
%     fh  high end of highest filter as a fraction of fs [default = 0.5]
%
%		w   any sensible combination of the following:
%
%               'R'  rectangular window in time domain
%				'N'	 Hanning window in time domain
%				'M'	 Hamming window in time domain (default)
%               'S'  Input is stft(time,freq), i.e. enframe and rfft have already been done
%
%               't'  triangular shaped filters in mel domain (default)
%               'f'  triangular filters in Hz domain using v_filtbankm instead of v_melbankm (better for large p)
%               'n'  hanning shaped filters in mel domain
%               'm'  hamming shaped filters in mel domain
%
%               'z'  highest and lowest filters taper down to zero (default)
%               'y'  lowest filter remains at 1 down to 0 frequency and
%			   	     highest filter remains at 1 up to nyquist freqency
%
%				'p'	 filters act in the power domain
%				'a'	 filters act in the absolute magnitude domain (default)
%
%               '0'  include 0'th order cepstral coefficient
%				'E'  include log energy
%				'd'	 include delta coefficients (dc/dt)
%				'D'	 include delta-delta coefficients (d^2c/dt^2)
%
%		       If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:	c     mel cepstrum output: one frame per row. Log energy, if requested, is the
%                 first element of each row followed by the delta and then the delta-delta
%                 coefficients.
%           tc    fractional time in samples at the centre of each frame
%                 with the first sample being 1.
%           y     log mel spectrogram output: one frame per row ('p' implies power domain)
%           mc    mel centre frequencies
%

% BUGS: (1) should have power limit as 1e-16 rather than 1e-6 (or possibly a better way of choosing this)
%           and put into VOICEBOX
%       (2) get v_rdct to change the data length (properly) instead of doing it explicitly (wrongly)
%       (3) need to sort out spectrogram plot x-axis when 'S' option is given

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: v_melcepst.m 10865 2018-09-21 17:22:45Z dmb $
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
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 fs=11025; end                       % default smple frequency = 11.025 kHz
if nargin<3 w='M'; end                          % default is Hamming time-domain window
if nargin<4 nc=12; end                          % default is 12 cepstral coefficient
if nargin<5 p=floor(3*log(fs)); end             % default mel-filterbank size (29 for fs=16kHz)
if nargin<6 n=pow2(floor(log2(0.03*fs))); end   % default frame length (30 ms)
if nargin<9
    fh=0.5;
    if nargin<8
        fl=0;
        if nargin<7
            inc=floor(n/2);
        end
    end
end
if isempty(w)
    w='M';
end
if any(w=='S')
    f=s.';                      % input s is stft(time,freq)
    tc=(1:size(f,2))';          % output frame number instead of centre-frequency
    if nargin<6                 % no frame length specified
        n=2*size(f,1)-2;            % assume original frame length was even
    end
else
    if any(w=='R')
        [z,tc]=v_enframe(s,n,inc);
    elseif any (w=='N')
        [z,tc]=v_enframe(s,0.5-0.5*cos(2*pi*(1:n)'/(n+1)),inc); % Hanning window
    else
        [z,tc]=v_enframe(s,0.54-0.46*cos(2*pi*(0:n-1)'/(n-1)),inc); % Hamming window
    end
    f=v_rfft(z.');              % f(freq,time) is STFT
end
if any(w=='f')
    if any(w=='y')
        [m,mc,a,b]=v_filtbankm(p,n,fs,fl,fh*fs,'mSxXH');
    else
        [m,mc,a,b]=v_filtbankm(p,n,fs,fl,fh*fs,'mSH');
    end
else
    wmb='tnmyz';                % list of input options relevant to v_melbankm
    [m,mc,a,b]=v_melbankm(p,n,fs,fl,fh,wmb(any(repmat(w',1,length(wmb))==repmat(wmb,length(w),1),1)));
end
pw=f(a:b,:).*conj(f(a:b,:));
pth=max(pw(:))*1E-20;
if any(w=='p')
    y=log(max(m*pw,pth));
else
    ath=sqrt(pth);
    y=log(max(m*abs(f(a:b,:)),ath));
end
c=v_rdct(y).';
nf=size(c,1);
nc=nc+1;
if p>nc
    c(:,nc+1:end)=[];
elseif p<nc
    c=[c zeros(nf,nc-p)];
end
if ~any(w=='0')
    c(:,1)=[];
    nc=nc-1;
end
if any(w=='E')
    c=[log(max(sum(pw),pth)).' c];
    nc=nc+1;
end

% calculate derivative

if any(w=='D')
    vf=(4:-1:-4)/60;
    af=(1:-1:-1)/2;
    ww=ones(5,1);
    cx=[c(ww,:); c; c(nf*ww,:)];
    vx=reshape(filter(vf,1,cx(:)),nf+10,nc);
    vx(1:8,:)=[];
    ax=reshape(filter(af,1,vx(:)),nf+2,nc);
    ax(1:2,:)=[];
    vx([1 nf+2],:)=[];
    if any(w=='d')
        c=[c vx ax];
    else
        c=[c ax];
    end
elseif any(w=='d')
    vf=(4:-1:-4)/60;
    ww=ones(4,1);
    cx=[c(ww,:); c; c(nf*ww,:)];
    vx=reshape(filter(vf,1,cx(:)),nf+8,nc);
    vx(1:8,:)=[];
    c=[c vx];
end

if nargout<1
    [nf,nc]=size(c);
    %    t=((0:nf-1)*inc+(n-1)/2)/fs;
    ci=(1:nc)-any(w=='0')-any(w=='E');
    imh = imagesc(tc/fs,ci,c.');
    axis('xy');
    xlabel('Time (s)');
    ylabel('Mel-cepstrum coefficient');
    map = (0:63)'/63;
    colormap([map map map]);
    colorbar;
elseif nargout>2
    y=y.'; % transpose y
end

