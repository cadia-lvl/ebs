function s=istfte(stft,meta,par)
% Epoch-based inverse stft
%
%  Inputs: stft(nframe,maxfft)     complex STFT coefficients (unused entries are normally NaN)
%          meta(nframe,2 to 6)     metadata: meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay]
%                                            default values:  dft-length=size(stft,2), offset=0, scale-factor=1, group-delay=0
%
% Outputs: s(ns,1)                  real-valued signal
%
% Transformations are applied in the order groupdelay, inverse dft, unpad, scale, offset
% Frames will be overlap-added to create the output, s, without any windowing.
%
% Bugs/Suggestions:
% (1) Add optional windowing before overlap-add if frames overlap.
%
[nframe,maxfft]=size(stft);             % number and maximum length of frames
nmeta=size(meta,2);                     % number of metadata parameters
if nmeta<6
    meta(:,6)=0;                        % set default group delay to zero
    if nmeta<5
        meta(:,5)=1;                    % set default scale factor to unity
        if nmeta<4
            meta(:,4)=0;                % set default offset to zero
            if nmeta<3
                meta(:,3)=maxfft;       % set default dft length to size(stft,2)
            end
        end
    end
end
dodct=nargin>2 && isfield(par,'groupdelay') && strcmpi(par.groupdelay,'dct');
framelens=min(meta(:,2),maxfft);        % datapoints in each frame in samples
s=zeros(meta(end,1:2)*[1;1]-1,1);       % space for reconstituted signal
frameend=meta(:,1)+framelens-1;         % calculate last datsapoint of each frame
% assume that frames might all be of different lengths (could probably speed up the constant-frame-length case)
for i=1:nframe
    nfft=meta(i,3);                     % DFT length of this frame
    if meta(i,6)~=0                     % apply group delay linear phase shift
        stft(i,1:nfft)=stft(i,1:nfft).*exp(-2i*pi/nfft*meta(i,6)*[0:ceil(nfft/2)-1 zeros(1,1-mod(nfft,2)) 1-ceil(nfft/2):-1]); % apply group delay (except to Nyquist frequency)
    end
    if dodct
        sfr=irdct(stft(i,1:nfft).');    % perform inverse DCT
    else
        sfr=ifft(stft(i,1:nfft));           % Perform inverse DFT
    end
    sfr=sfr(1:framelens(i));            % truncate to remove any padding
    sfr=sfr*meta(i,5)+meta(i,4);        % scale the data and add offset
    s(meta(i,1):frameend(i))=s(meta(i,1):frameend(i))+real(sfr(:)); % overlap add into the output waveform
end
if ~nargout
    fredge=unique([meta(:,1)-0.5; meta(:,1:2)*[1;1]-0.5]); % frame boundaries
    plot(s);
    hold on;
    if ~all(meta(:,4)==0)               % if there are any non-zero offsets
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[1;1]*meta(:,4)',':k');
    end
    if ~all(meta(:,5)==1)               % if there are any non-unity scale factors
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[1;1]*meta(:,5)','--r');
        plot([meta(:,1) meta(:,1:2)*[1;1]-1]',[-1;-1]*meta(:,5)','--r');
    end
    v_axisenlarge([-1 -1.05]);
    ylim=get(gca,'ylim');
    plot([1;1]*fredge',ylim','k:');     % mark frame boundaries
    if ~all(meta(:,6)==0)               % if there are any non-zero group delays
        plot(meta(:,[1 6])*[1;1]-1,ylim*[0.97;0.03],'^r',meta(:,[1 6])*[1;1]-1,ylim*[0.03;0.97],'vr');
    end
    hold off
    xlabel('Sample Number');
    ylabel('Signal');
end