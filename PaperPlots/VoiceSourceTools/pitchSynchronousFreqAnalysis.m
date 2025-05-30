function [spseg] = pitchSynchronousFreqAnalysis(sp, fs, gci)

%pitchSynchronousFreqAnalysis 
%
%  Inputs:  sp       is the speech signal 
%           fs       is the sampling frequency in Hz
%           gci      is the glottal closure instants indices (samples) of
%                    length Ng
%
% Outputs: 
% 

nfft = 2^12;

for ig=2:length(gci)-1
    %Create a window around the GCI.  The current implementation uses the
    %average of the before and after pitch periods as a window size, but it
    %could be set arbitarily
    %wz=(gci(ig+1)-gci(ig-1));
    %wz=600;
    %nn=(gci(ig)-ceil(wz/2)):gci(ig)+ceil(wz/2);%+round(0.5*wz*rand);
    nn=gci(ig-1):gci(ig);
    wz=length(nn);
    
    f0=fs/(gci(ig)-gci(ig-1));
    
      
    % If the indices exceed the signal boundaries then cut 
    nn=nn(nn>0); 
    nn=nn(nn<length(sp));

    spseg=sp(nn);
    Sf=fft(hamming(wz).*spseg,nfft);
    ff=fs*(0:nfft-1)/nfft;
    Sf=Sf(1:nfft/2);
    ff=ff(1:nfft/2);
    
    figure(1); clf; set(gcf,'Position', [12   385   560   420])
    plot(spseg);
    title(['Instantaneous f0 : ' num2str(f0) ' Hz'])
    figure(2); clf; set(gcf,'Position', [788    68   560   737])
    ax(1)=subplot(211);    
    plot(ff,10*log10(abs(Sf)));
    %xlim([0 5000])
    ax(2)=subplot(212);
    plot(ff,unwrap(angle(Sf))/(2*pi*fs));
    ylabel('Phase [s]')
    %xlim([0 5000])
    linkaxes(ax,'x');
    
    ufe=2;
    

    %Spectral parameters
    %[pxx1,f] = pwelch(useg-mean(useg),hann(T),0,nfft,fs,'onesided');
    %PSD1=10*log10(pxx1/max(pxx1));
    
%     figure(1); clf;
%     plot(spseg);
%     pause
    
   
end;
