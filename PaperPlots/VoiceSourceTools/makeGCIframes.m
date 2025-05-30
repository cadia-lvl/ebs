function [frameStart, vu, gci, gciStruct] = makeGCIframes(sp, fs, timeMarkers, gci, UVmethod)

%MAKEGCIFRAMES extract GCIs and give framestarts from UV segmented speech [frameStart, gci, gciStruct] = (sp, fs, timeMarkers)
%
%  Inputs:  sp     is the input signal
%           fs     is the sampling frequency
%           timeMarkers  specifies the voiced/unvoiced segmentation
%
% Outputs:  frameStart  The start of analysis frames (regular in U, GCIs in V)
%           vu          Voiced=1, Unvoiced=0 for each frame
%           gci         Glottal closing instants
%           gciStruct   
%           

%Parameters 
Nlen=length(sp);
minsz=fs*1e-3;
maxsz=fs*20e-3;
winsz=fs*30e-3;         % Frame size in unvoiced segments
smallTh=fs*30e-3;       % If a voiced segment is too small (unvoiced segments need to be larger than winsz)

% Dilate and make sure the segments do not index out of the signal
nExt= 0;%round(10e-3*fs);  % Extend (dilate) the voiced segment by 10 ms either side

iT=max(round(timeMarkers(:,1)*fs)-nExt,1);
eT=min(round(timeMarkers(:,2)*fs)+nExt,Nlen);

if nargin<4
    gci=[];
end

if nargin < 5
    UVmethod='GCI';
end

frameStart=1;
vu=0;
   
for ii=1:length(iT)
    % The first voiced segment
    segStart = iT(ii);
    segEnd   = eT(ii);
    
    if (segEnd-segStart)< smallTh   % Do not detect GCIs in small segments
        %gciAbs=segStart:winsz:segEnd; %
        gciAbs=segEnd;
    else
        if nargin<4
            % GCI detection : would be nice to dynamically control
            % voicebox('dy_fxmax',500); and voicebox('dy_fxmin', 50);
            gciSeg = dypsagoi(sp(segStart:segEnd),fs);
        else
            gciSeg=gci((segStart <= gci & gci <= segEnd))-segStart+1;
        end
        gciAbs = gciSeg+segStart-1;
    end
    % Checking GCIs - there should be a switch in dypsagoi (max f0)
    % numSmall(ii)=sum(1000*diff(gciSeg)/fs<2);
    
    if nargout>3
        %gci(end+1:end+length(gciSeg)) = gciAbs;
        gciStruct(ii).gciSeg = gciSeg;
        gciStruct(ii).segStart = segStart;
        gciStruct(ii).segEnd = segEnd;
    end
    
    % Add unvoiced frames in the segment preceding the current voiced
    % part
    switch UVmethod
        case 'Fixed'
            unVoicedFrames=(frameStart(end)+winsz):winsz:(gciAbs(1)-round(winsz/2));
        case 'GCI'
            unVoicedFrames=gci((frameStart(end)+1 <= gci & gci <= segStart));
            if length(unVoicedFrames)<3
                unVoicedFrames=(frameStart(end)+minsz):minsz:(gciAbs(1)-round(minsz/2));
            end
            if (unVoicedFrames(1) - frameStart(end)) < minsz
                unVoicedFrames(1)=[];
            end
            if (gciAbs(1) - unVoicedFrames(end)) < minsz
                unVoicedFrames(end)=[];
            end
            idxsmall=diff(unVoicedFrames)<minsz;
            if any(idxsmall)
                unVoicedFrames(idxsmall) = [];
            end
            idxbig=find(diff(unVoicedFrames)>maxsz);
            for ix=1:length(idxbig)
                mgci=round((unVoicedFrames(idxbig(ix))+unVoicedFrames(idxbig(ix)+1))/2);
                unVoicedFrames = [unVoicedFrames(1:idxbig(ix)) mgci unVoicedFrames((idxbig(ix)+1):end)];
                idxbig(ix+1:end)=idxbig(ix+1:end)+1;
            end
                
                
    end
    frameStart=[frameStart unVoicedFrames gciAbs];
    vu=[vu zeros(size(unVoicedFrames)) ones(size(gciAbs))];
    
    
end
%Last unvoiced segment
    switch UVmethod
        case 'Fixed'
            unVoicedFrames=(frameStart(end)+winsz):winsz:(Nlen-round(winsz/2));
        case 'GCI'
            unVoicedFrames=gci((frameStart(end)+1 <= gci & gci <= Nlen));
    end
    
frameStart=[frameStart unVoicedFrames];
vu=[vu zeros(size(unVoicedFrames))];

if nargout==0
    figure(22); clf;
    plot(1000*diff(frameStart)/fs); hold on;
    xlabel('FrameIx')
    ylabel('Time [ms]')
    title('Length of windows (diff(frameStart))')
    
    figure(19); clf;
    stem(frameStart(vu==1), 0.3*ones(size(frameStart(vu==1))))
    hold on
    stem(frameStart(vu==0), 0.3*ones(size(frameStart(vu==0))))
    plot(sp)
end

%figure(23); clf;
%plot(numSmall);



    
    
