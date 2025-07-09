function [framek,gci]=gs_frames(s,fs,p)
% Divide speech signal into pitch-synchronous frames
%
%  Inputs: s(nsamp,1)       Speech signal
%          fs               Sample Frequency
%          p                optional structure giving procesing options:
%                               p.pitchlim       vector with [min target max] pitch (Hz)
%                               p.gcifrac        position of GCI in analysis frame [0.3]
%                               p.GCImethod      either 'YAGA' (default) or 'SEDREAMS'
%
% Outputs: framek(1,nframe) Vector containng the index of the last sample in each frame
%          gci(1,ngci)      Glottal Closure Instants (seconds, origin @ s(1)) (dummy GCIs are inserted if necessary to make spacing <= 1/min-pitch)
%
% If the GCI option p.GCImethod='SEDREAMS' is selected, this routine uses pitch_srh and
% gci_sedreams from covarep. Note that after installing covarep, you must remove its version
% of voicebox from the MATLAB path (because it is out-of-date).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   algorithm Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent q0
if isempty(q0)
    q0.pitchlim=[40 50 400];                               % {min target max} pitch (Hz)
    q0.gcifrac=0.3;                                        % position of GCI in analysis frame
    q0.GCImethod='YAGA';
end
if nargin<3
    q=q0;
else
    q=v_paramsetch(q0,p);             % update parameters from p input
end
%
ns=length(s);                                           % length of speech signal
switch q.GCImethod                                      % find gci's
    case 'YAGA'
        [gci,goi] = dypsagoi(s,fs);
        gci=gci/fs;
    case 'SEDREAMS'
        [fx,vad,srh,fxt] = pitch_srh(s,fs,q.pitchlim(1),q.pitchlim(3));     % find pitch track
        fxmed=median(fx);                               % median pitch
        [gci, mbs, res] = gci_sedreams(s, fs, fxmed, 1);
    otherwise
        error(['GCI method ' q.GCImethod ' is not recognised.']);
end;
gcik=round(gci*fs+1);                                   % convert to sample numbers
gcikorig=gcik;                                          % save for debugging
periodlimk=round(fs./q.pitchlim);                       % {max targt min} pitch periods (samples)
%
% %%%%%%%%%%%%% Note: this could probably all be speeded up
% %%%%%%%%%%%%% might need to delete first few or last few gcis if they are too close to the ends
%
% delete very short pitch periods
%
mk=(gcik(2:end)-gcik(1:end-1))<periodlimk(3);           % find pitch periods that are too short
gcid=[]; % space for deleted gcis
while any(mk)
    mkd=[false mk] & [false true ~mk(1:end-1)];         % identify the first of a string of short periods
    gcid=[gcid gcik(mkd)];                              % append to list of deleted GCIs
    gcik(mkd)=[];                                       % delete last gci in the first of each string of short periods
    mk=(gcik(2:end)-gcik(1:end-1))<periodlimk(3);       % find pitch periods that are too short
end
%
% insert dummy gcis into very long pitch periods
%
mk=(gcik(2:end)-gcik(1:end-1))>periodlimk(1);
bigpp=find(mk);
gaps=gcik(bigpp+1)-gcik(bigpp);
nins=ceil(gaps/periodlimk(2))-1;                        % number of gcis to insert
incs=gaps./(nins+1);                                    % desired increment (fractional samples)
ninsc=[0 cumsum(nins)];
gcix=zeros(1,ninsc(end));                               % space for added gcis
for i=1:length(nins)
    gcix(ninsc(i)+1:ninsc(i+1))=gcik(bigpp(i))+round((1:nins(i))*incs(i));
end
%
% Note: dft frame is GCI - round(q.gcifrac*len) + (1:len), gci(i)-round(q.gcifrac*(gci(i)-gci(i-1)))+1:gci(i+1)-round(q.gcifrac*(gci(i+1)-gci(i)))
%
% calculate gcis to insert at the start
nins=ceil(gcik(1)/periodlimk(2)-q.gcifrac+1);                     % number of gcis to insert at start including a negative one
incs=gcik(1)/(nins+q.gcifrac-1);                                  % desired increment (fractional samples)
gciy=gcik(1)-round((1:nins)*incs);                              % inserted GCIs at the start
% now do gcis at the end
nins=ceil((ns-gcik(end))/periodlimk(2)+q.gcifrac);                % number of gcis to insert at end including one beyond ns
incs=(ns-gcik(end))/(nins-q.gcifrac);
gciz=gcik(end)+round((1:nins)*incs);
% add the inserted gcis
gcia=[gcix gciy gciz];
gcim=sort([gcik gcia]);
framek=gcim(3:end)-round(q.gcifrac*(gcim(3:end)-gcim(2:end-1)));  % frame end samples