function [framek,gci]=gs_frames(s,fs,pardb)
% Divide speech signal into pitch-synchronous frames
%
%  Inputs:      s       Speech signal
%              fs       Sample Frequency
%               pardb (optional)    
%                         .pitchlim   {min target max} pitch (Hz)
%                         .gcifrac    position of GCI in analysis frame
%                         .GCImethod  either 'YAGA' (default) or SEDREAMS
%
% Outputs: framek   Vector containng the index of the last sample in each frame
%           gci     Glottal Closure Instants
%
% This routine uses pitch_srh and gci_sedreams from covarep. Note that after installing covarep,
% you must remove its version of voicebox from the MATLAB path (because it is out-of-date).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   algorithm Parameters
%   (could use an input parameter structure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    pitchlim=[40 50 400];                               % {min target max} pitch (Hz)
    gcifrac=0.3;                                        % position of GCI in analysis frame
    GCImethod='YAGA';
else
    pitchlim=pardb.pitchlim;                               % {min targt max} pitch (Hz)
    gcifrac=pardb.gcifrac;                                 % position of GCI in analysis frame
    GCImethod=pardb.GCImethod;
end
%
ns=length(s); % length of speech signal
switch GCImethod % find gci's
    case 'YAGA'
        [gci,goi] = dypsagoi(s,fs);
        gci=gci/fs;
    case 'SEDREAMS'
        [fx,vad,srh,fxt] = pitch_srh(s,fs,pitchlim(1),pitchlim(3));     % find pitch track
        fxmed=median(fx);                                               % median pitch
        [gci, mbs, res] = gci_sedreams(s, fs, fxmed, 1);
    otherwise
        error(['GCI method ' GCImethod ' is not recognised.']);
end;
gcik=round(gci*fs+1);                                           % convert to sample numbers
gcikorig=gcik; % save for debugging
periodlimk=round(fs./pitchlim);                         % {max targt min} pitch periods (samples)
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
    gcid=[gcid gcik(mkd)]; % append to list of deleted GCIs
    gcik(mkd)=[];    % delete last gci in the first of each string of short periods
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
% Note: dft frame is GCI - round(gcifrac*len) + (1:len), gci(i)-round(gcifrac*(gci(i)-gci(i-1)))+1:gci(i+1)-round(gcifrac*(gci(i+1)-gci(i)))
%
% calculate gcis to insert at the start
nins=ceil(gcik(1)/periodlimk(2)-gcifrac+1);                     % number of gcis to insert at start including a negative one
incs=gcik(1)/(nins+gcifrac-1);                                  % desired increment (fractional samples)
gciy=gcik(1)-round((1:nins)*incs);                              % inserted GCIs at the start
% now do gcis at the end
nins=ceil((ns-gcik(end))/periodlimk(2)+gcifrac);                % number of gcis to insert at end including one beyond ns
incs=(ns-gcik(end))/(nins-gcifrac);
gciz=gcik(end)+round((1:nins)*incs);
% add the inserted gcis
gcia=[gcix gciy gciz];
gcim=sort([gcik gcia]);
framek=gcim(3:end)-round(gcifrac*(gcim(3:end)-gcim(2:end-1)));  % frame end samples