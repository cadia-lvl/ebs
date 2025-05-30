function ax=plotFourierCoefficients(fh,sSeg,ffrange)

if nargin<3
    ffrange=[0 20000];
end

fh=figure(fh); clf;
ax(1)=subplot(211);
plot(sSeg.ff,10*log10(abs(sSeg.a)));
ylabel('Log-magnitude');
ax(2)=subplot(212);
plot(sSeg.ff,1000*sSeg.N*unwrap(angle(sSeg.a))/(sSeg.fs*2*pi));
ylabel('Phase shift [ms]');
xlabel('Frequency [Hz]');

linkaxes(ax,'x')
xlim(ffrange);