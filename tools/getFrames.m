function frames=getFrames(s,fs,par,frameMethod);

% getFrames Calculates frame boundaries for speech spectrogram calculation
%
%  The function selects a method to extract frame boundaries from the
%  parameter par.frameMethod{ik}
%
%
%  Inputs:  s    Speech signal vector
%           fs   Sampling frequency of the speech 
%           par      Parameter structure containing:
%              .frameMethod
%              .nfft
%              .maxoff
%              .periodlimk
%
% Outputs:  frames  An [Nf x 2] dimensional vector containing the first and
%                   the last sample value of each of Ns frames for the
%                   signal

ns=length(s);
switch frameMethod
    case 'Fixed'
        frames=par.nfft:par.nfft:ns;
    case 'FixedNE'
        frames=gs_frames(s,fs,par);    %Just to figure out the number of frames of Epoch
        nfr=length(frames);
        nfft=round(ns/nfr);
        frames=nfft:nfft:ns;
    case 'Epoch'
        frames=gs_frames(s,fs,par);
    case 'EpochAdj'
        framek=gs_frames(s,fs,par);
        frames=gs_frameadj(framek,s,fs,par.maxoff);
        framekk=[[1 frames(1:end-1)+1];frames];             %  start and end samples for each frame
        if any(framekk(:)<1) || any(framekk(:)>ns)
            error('frame boundaries do not lie within 1:ns');
        end
        framelen=framekk(2,:)-framekk(1,:)+1;
        if (any(framelen>par.periodlimk(1)))
            warning('Some frames are too long')
        end
    otherwise
        error(['Frame method " ' frameMethod ' "is not defined']);

end


