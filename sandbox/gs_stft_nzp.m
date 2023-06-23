function [stft,framelen]=gs_stft_nzp(s,framex,nfft)

nframe=length(framex); % number of frames
framekk=[[1 framex(1:end-1)+1];framex]; %  start and end samples for each frame
framelen=framekk(2,:)-framekk(1,:)+1;
stft=zeros(nfft,nframe);

for i=1:nframe
    stft(1:framelen(i),i)=fft(s(framekk(1,i):framekk(2,i)));
end
