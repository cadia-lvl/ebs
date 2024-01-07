function [x,y,b]=fftgdopterr(g,v,w,a)
% calculate phase error for optimization to find the best group delay
% needed for the par.groupdelay='fmbd' option to stfte()
%
% Usage:    dft=fft(...);                                         % complex-valued dft vector
%           [x,y,b]=fftgdopterr(dft);                             % initial call to compute quantities independent of group delay
%           gopt = fminbnd(@(g) fftgdopterr(g,x,y,b),0,nfft-1);   % find optimal group delay
%
%  Inputs: g    If nargin>1: group delay in samples; normally in the range 0 to nfft-1 where nfft is the dft length
%               If nargin=1: complex-valued dft vector
%          v    the x output from the intial call to fftgdopterr(dft)
%          w    the y output from the intial call to fftgdopterr(dft)
%          a    the b output from the intial call to fftgdopterr(dft)
%
% Outputs: x    If nargin>1: the energy-weighted phase error
%               If nargin=1: complex-valued dft vector; positive frequencies only excluding DC and Nyquist
%          y    If nargin=1: scaled absoluted value of the dft vector; positive frequencies only excluding DC and Nyquist
%          b    If nargin=1: scalar constant: 2i*pi/nfft where nfft is the dft length
%
% This routine is used by the fminbnd() optimizer to find the group delay that best matches the phase spectrum of a DFT.
% The metric that is minimized is an energy-weighted average of 1-cos(phi) where phi is an element of the phase spectrum.
% This metric lies in the range 0 to 2.
%
% energy-weighted error = sum(abs(x)^2*(1-cos(angle(x)))) / sum(abs(x)^2)
%                       = 1 - sum(abs(x)^2*cos(angle(x)))) / sum(abs(x)^2)
%                       = 1 - sum(abs(x)*real(x)) / sum(abs(x)^2)
%
if nargin>1                 % called from the optimizer so need to calculate the error
    x=1-w'*real(v.*cumprod(repmat(exp(a*g),numel(v),1))); % calculate error
else                        % initial call to precalculate fixed quantities (indicated by only one input argument)
    g=g(:);                 % force to be a column vector
    nfft=numel(g);          % dft length
    nfftq=ceil(nfft/2);     % number of positive frequencies excluding Nyquist
    x=g(2:nfftq);           % positive frequencies excluding DC and Nyquist
    y=abs(x);               % also need absolute value
    y=y/(y'*y);             % divide by total energy
    b=2i*pi/nfft;           % scaling constant converts group delay from samples to phase increment
end
