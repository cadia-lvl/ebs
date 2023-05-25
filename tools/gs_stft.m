function [stft,fax]=gs_stft(s,framex,nfft)

% GS_STFT Calculates the Short Time Fourier Transform of s based on
% preefined frames
%
%  Inputs:     s        Speech signal
%              framex   Frame positions
%              nfft     Number of Fourier coefficients to calculate for each frame   
%
% Outputs:     stft     The Short Time Fourier Transform
%              fax      Frequency (in Hz) of each frequency bin
%

nframe=length(framex); % number of frames
framekk=[[1 framex(1:end-1)+1];framex];                         %  start and end samples for each frame
framelen=framekk(2,:)-framekk(1,:)+1;
sfr=zeros(nfft,nframe);
frame0=1;
for i=1:nframe
    sfr(1:framelen(i),i)=s(framekk(1,i):framekk(2,i));
    sfr(framelen(i)+1:nfft,i)=[0.5 0.5]*sfr([1 framelen(i)],i); % pad with average of end points
end
stft=fft(sfr);                                                  % calculate dft
nfftp=1+floor(nfft/2);                                          % number of positive frequencies
fax=(0:nfftp-1)/nfft;