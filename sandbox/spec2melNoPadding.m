function stftm=spec2melNoPadding(stft,fs,numMelFilters,framelen)

%  Inputs: sdft(nfft,nframes)    complex STFT coefficients
%          fs                   sample frequency                                
%          numMelFilters        number of mel filters
%          framelen(1,nframes)  Frame lengths of each frame.  The mel-filters assume that
%                               the spectrum is contained in the first "framelen(1,i)" samples. 
%
% Outputs: stftr(nfft,nframe)           Reconstructed complex STFT
%          stftm(numMelFilters,nframe)  Real-valued Mel filterbank outputs

nframes=size(stft,2);

stftm=zeros(numMelFilters,nframes);
for ii=1:nframes
    nfftp=1+floor(framelen(ii)/2);  
    mbm=v_filtbankm(numMelFilters,framelen(ii),fs,0,fs/2,'m');
    stftm(:,ii)=mbm*abs(stft(1:nfftp,ii).^2);
end
