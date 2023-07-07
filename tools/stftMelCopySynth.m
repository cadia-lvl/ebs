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

for ik=1:length(par.config)
    frameMethod=par.config{ik}.frameMethod;
    preserveDC=par.config{ik}.preserveDC;
    nzp=par.config{ik}.nzp;

    % 1) Calculate frames for method ik
    frames=getFrames(s,fs,par,frameMethod);

    % 2) Calculate STFT based on frames
    if nzp
        [stft,framelen]=gs_stft_nzp(s,frames,par.nfft);    
    else
        stft=gs_stft(s,frames,par.nfft);
    end

    for is=1:length(par.numMels)
        % 3) Convert spectrograms to Mel energy - grams
        % 4) Convert Mel Energy-grams to spectrograms
        % 5) Invert spectrograms to time domain signal
        if nzp
           stftm=spec2mel_nzp(stft,fs,par.numMels(is),framelen,preserveDC);
           stftr=mel2spec_nzp(stftm,fs,par.nfft,framelen,preserveDC,angle(stft));
           vsr=gs_istft_nzp(stftr,frames);
        else
           stftm=spec2mel(stft,fs,par.numMels(is),preserveDC);
           stftr=mel2spec(stftm,fs,par.nfft,preserveDC,angle(stft));
           vsr=gs_istft(stftr,frames);
        end

        % 6) Calculate SNR and maybe PESQ
        snrMat(is,ik)=calculateSNR(s,vsr);

        if par.doPESQ==1
            pesqMat(is,ik)=computePESQ(s,vsr,par);
        end

        % Calculate SNR for each phone segment
        nsu=length(vsr);
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


