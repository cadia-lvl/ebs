function sSeg=gCycleSegments(CIx,sp,fs,gci,PhaseShift,wname,Nfft)


if nargin < 5
    PhaseShift=0;
    if nargin < 6
        wname=@rectwin;
    end
end

tt=(0:length(sp)-1)/fs;

nn=(gci(CIx(1)):(gci(CIx(end)+1)-1))-PhaseShift;
sp=sp(nn); 
ttSeg=tt(nn);
N=length(sp);
if nargin < 7
    Nfft=N;
end
a=fft(sp.*window(wname,N),Nfft)/N;
ff=(fs/Nfft)*(0:(Nfft-1));

sSeg=Var_Names(CIx,nn,sp,ttSeg,N,a,ff,PhaseShift,wname,fs,Nfft);


function out = Var_Names(varargin)
for n = 1:nargin
    out.(inputname(n)) = varargin{n};
end