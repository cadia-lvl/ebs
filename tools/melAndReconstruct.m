function [sdftr, sdftm]=melAndReconstruct(sdft,fs,nfft,nfftp,numMelFilters,fbankmethod,preserveDC)
%  Inputs: sdft(nfft,nframe)    complex STFT coefficients
%          fs                   sample frequency
%          nfft                 FFT length
%          nfftp                Num of positive frequencies = 1+floor(nfft/2)
%          numMelFilters        number of mel filters
%          method               Can be 'melbankm' or 'filtbankm'
%          preserveDC           If set to 1, the DC value of sdft is copied
%                               into sdftr.  Otherwise, estimated DC values are used 
%
% Outputs: sdftr(nfft,nframe)  Reconstructed complex STFT
%          sdftm(numMelFilters,nframe) Real-valued Mel filterbank outputs

if nargin < 7
    preserveDC = 0;
    if nargin < 6
        fbankmethod='filtbankm';
    end;
end;

switch fbankmethod
    case 'filtbankm'
        mbm=v_filtbankm(abs(numMelFilters),nfft,fs,0,fs/2,'m');
    case 'melbankm'
        mbm=v_melbankm(abs(numMelFilters),nfft,fs);
    otherwise
        error(['fbankmethod ' fbankmethod ' not defined'])
end
if preserveDC
    mbm=vertcat(sparse(1,1,1,1,nfftp),mbm); % preappend an extra row to preserve the DC value
end
imbm=pinv(full(mbm));
sdftm=mbm*abs(sdft(1:nfftp,:).^2);
sdftr=sqrt(max(imbm*sdftm,0)).*exp(1i*angle(sdft(1:nfftp,:)));
sdftr=[sdftr; conj(sdftr(end-1:-1:2,:))];

