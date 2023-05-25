function out = stftMelCopySynth(spfile,par)

% stftMelCopySynth performs spectrogram Mel copy-synthesis
%
%  The function takes a (TIMIT) speech file name, computes spectrograms as
%  described in par.frameMethod, parameterises it using a number of
%  Mel-filters (defined in par.numMels, possibly many) and returs
%  signal-to-noise ratios (and PESQ scores if required) for each method,
%  choice of number of Mel-filters and for each phoneme (TIMIT)
%
%  Inputs:  spfile   Speech file name (wav spherical format)
%           par      Parameter structure containing:
%              .numMels
%              .frameMethod
%              .doPESQ
%              .nfft
%
% Outputs:  out     An output structure containing:
%              .snrMat
%              .snrMatPhn
%              .pesqMat  (optional) 
%              .phn
%              .spfile

[s,fs]=v_readsph(spfile,'wt');
phn=extractPhn([spfile(1:end-3) 'PHN']);

snrMat=zeros(length(par.numMels),length(par.frameMethod));

if par.doPESQ==1
    pesqMat=zeros(length(par.numMels),length(par.frameMethod));
end

snrMatPhn=cell(length(par.numMels),length(par.frameMethod),length(phn.in));

for ik=1:length(par.frameMethod)
    % Calculate frames for method ik
    frames=getFrames(s,fs,par,ik);
    
    % Calculate STFT based on frames
    [sdft,fax]=gs_stft(s,frames,par.nfft);
    fax=fax*fs;                                         % convert frequency axis to Hz
    nfftp=length(fax);

    for is=1:length(par.numMels)
        % Calculate mel spectrum and reconstruct
        sdftr=melAndReconstruct(sdft,fs,par.nfft,nfftp,par.numMels(is));

        %Do the inverse FFT
        vsr=gs_istft(sdftr,frames);

        nsu=length(vsr);                                 % lengthof vu is slightly less due to framing
        ssq=sum(s(1:nsu).^2);
        snrMat(is,ik)=-db(sum((vsr(1:nsu)-s(1:nsu)).^2)/ssq);

        if par.doPESQ==1
            pesqMat(is,ik)=computePESQ(s,vsr,par);
        end

        for ij=1:length(phn.tstart)
            segIx=max([1,phn.tstart(ij)]):min([nsu,phn.tend(ij)]);
            ssqSeg=sum(s(segIx).^2);
            snrMatPhn{is,ik,ij}=-db(sum((vsr(segIx)-s(segIx)).^2)/ssqSeg);
        end

    end
end
out.snrMat=snrMat;
if par.doPESQ==1
    out.pesqMat=pesqMat;
end
out.snrMatPhn=snrMatPhn;
out.phn=phn;
out.spfile=spfile;


