function [sdftr, sdftm]=melAndReconstruct(sdft,fs,nfft,nfftp,numMelFilters,method)
%  Inputs: sdft(nfft,nframe)    complex STFT coefficients
%          fs                   sample frequency
%          nfft                 FFT length
%          nfftp                Num of positive frequencies = 1+floor(nfft/2)
%          numMelFilters        +- number of mel filters (DC wil be preserved if negative)
%          method               Can be 'melbankm' or 'filtbankm'
%
% Outputs: sdftr(nfft,nframe)  Reconstructed complex STFT
%          sdftm(numMelFilters,nframe) Real-valued Mel fiolterbank outputs


if nargin>5 && strcmpi(method,'filtbankm')
    mbm=v_filtbankm(abs(numMelFilters),nfft,fs,0,fs/2,'m');
else
    mbm=v_melbankm(abs(numMelFilters),nfft,fs);
end
imbm=pinv(full(mbm));
sdftm=mbm*abs(sdft(1:nfftp,:).^2);
sdftr=sqrt(max(imbm*sdftm,0)).*exp(1i*angle(sdft(1:nfftp,:)));
sdftr=[sdftr; conj(sdftr(end-1:-1:2,:))];