function c=mel2melcep(y,w,nc)
% Inputs:
%     y(nt,np)	  real-valued mel spectrogram (nf frames, np filters)
%     w   mode string (see below)
%     nc  number of cepstral coefficients excluding 0'th coefficient [default 12]
%
%		w   any sensible combination of the following:
%
%               '0'  include 0'th order cepstral coefficient
%				['E'  include log energy] not yet implemented
%				'd'	 include delta coefficients (dc/dt)
%				'D'	 include delta-delta coefficients (d^2c/dt^2)
%
%		       If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:	c     mel cepstrum output: one frame per row. Log energy, if requested, is the
%                 first element of each row followed by the delta and then the delta-delta
%                 coefficients.
%
% Bugs/Suggestions:
%   (1) Needs to incorporate 'E' option and deal correctly with the
%       situation where v_filtbankm included the 'q','z' or 'Z' options
%
if nargin<3
    nc=12;                      % default to 12 cepstral coefficients
    if nargin<2
        w='';                   % dummy mode options
    end
end
pth=max(y(:))*1E-20; % low threshold to avoid infinite logs
c=v_rdct(log(max(y,pth)).').';
[nf,p]=size(c); % number of frames and filters
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
% if any(w=='E') % include log energy
%    c=[log(max(sum(pw),pth)).' c];
%    nc=nc+1;
% end

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