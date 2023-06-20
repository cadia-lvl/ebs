function v=gs_istft_nzp(stft,framex)
% No zero padding

framekk=[[1 framex(1:end-1)+1];framex];        %  start and end samples for each frame
framelen=filter([1,-1],1,framekk(2,:)); % frame lengths
v=zeros(sum(framelen),1); % space for reconstituted speech
ix=1; % frame start sample
for ii=1:length(framelen) % loop for each frame
    ilen=framelen(ii); % length of this frame
    jx=ix+ilen-1; % last sample from this frame
    v(ix:jx)=real((ifft(stft(1:ilen,ii)))); % do inverse transform and force real
    ix=jx+1;
end         