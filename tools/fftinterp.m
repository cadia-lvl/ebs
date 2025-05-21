function t=fftinterp(nin,nout)
% linearly interpolate dft spectrum onto new spacing
% output is a sparse interpolation matrix, t(nout,nin)
% this routine is inefficient because the interpolation coefficients for negative frequencies will always be a mirror image of those for positive frequencies
fin=(0:nin);
fout=(0:nout-1)*(nin/nout);
[foi,ix,jx]=v_sort([fin fout]);
iix=jx(nin+2:end)-(1:nout);             % index to next lower input frequency
frac=(fout-fin(iix));               % fraction of input hop above the lower frequency
t=sparse(repmat(1:nout,1,2),[iix mod(iix,nin)+1],[1-frac frac],nout,nin);