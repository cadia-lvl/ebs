function stftm=spec2mel(stft,fs,numMelFilters,preserveDC)

%  Inputs: sdft(nfft,nframes)    complex STFT coefficients
%          fs                   sample frequency                                
%          numMelFilters        number of mel filters
%          framelen(1,nframes)  Frame lengths of each frame.  The mel-filters assume that
%                               the spectrum is contained in the first "framelen(1,i)" samples. 
%
% Outputs: stftr(nfft,nframe)           Reconstructed complex STFT
%          stftm(numMelFilters,nframe)  Real-valued Mel filterbank outputs

if nargin < 4
    preserveDC = 0;
end;

nfft=size(stft,1);
nfftp=1+floor(nfft/2);
mbm=v_filtbankm(abs(numMelFilters),nfft,fs,0,fs/2,'m');
if preserveDC
    mbm=vertcat(sparse(1,1,1,1,nfftp),mbm); % preappend an extra row to preserve the DC value
end
stftm=mbm*abs(stft(1:nfftp,:).^2);
