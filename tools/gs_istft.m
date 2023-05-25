function v=gs_istft(stft,framex)

% GS_ISTFT Calculates the inverse of the Short Time Fourier Transform of s 
% based on preefined frames
%
%  Inputs:     stft     Short time Fourier transform
%              framex   Frame positions
%
% Outputs:     v        Estimated time signal
%

framekk=[[1 framex(1:end-1)+1];framex];        %  start and end samples for each frame
vfr=real(ifft(stft)); % do inverse transform and force real
vlen=filter([1,-1],1,framekk(2,:)); % frame lengths
v=zeros(sum(vlen),1); % space for reconstituted speech
ix=1; % frame start sample
for i=1:length(vlen) % loop for each frame
    ilen=vlen(i); % length of this frame
    jx=ix+ilen-1; % last sample from this frame
    v(ix:jx)=vfr(1:ilen,i);
    ix=jx+1;
end                