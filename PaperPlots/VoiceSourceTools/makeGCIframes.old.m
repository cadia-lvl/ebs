function [frameStart, gci, gciStruct] = makeGCIframes(sp, fs, timeMarkers)

%Parameters 
winsz=fs*30e-3;         % Frame size in unvoiced segments
smallTh=fs*50e-3;       % If a voiced segment is too small (unvoiced segments need to be larger than winsz)
nExt= 0;   %round(20e-3*fs);  % Extend (dilate) the voiced segment by 20 ms either side
Nlen=length(sp);

% Dilate and make sure the segments do not index out of the signal
iT=max(round(timeMarkers(:,1)*fs)-nExt,1);
eT=min(round(timeMarkers(:,2)*fs)+nExt,Nlen);

gci=[];

frameStart=1; % (Gets overwritten if there is an unvoiced segment at start)

% Do the first unvoiced segment (and skip if there isn't one)
if iT(1) > 1
    segStart=1;
    segEnd=iT(1)-1;
    frameStart=segStart:winsz:segEnd;
end
    
    
for ii=1:length(iT)
    segStart = iT(ii);
    segEnd   = eT(ii);
    
    if (segEnd-segStart)< smallTh   % Do not detect GCIs in small segments
        frameStart=[frameStart segStart]; % 
    else
        gciSeg = dypsagoi(sp(segStart:segEnd),fs);
        gciAbs = gciSeg+segStart-1;
        
        numSmall(ii)=sum(1000*diff(gciSeg)/fs<2);
        
        if nargout>1
            gci(end+1:end+length(gciSeg)) = gciAbs;
            gciStruct(ii).gciSeg = gciSeg;
            gciStruct(ii).segStart = segStart;
            gciStruct(ii).segEnd = segEnd;
        end
        
        % If the gap between last unvoiced frame start and the first GCI in
        % the new segment is too large, then split the frame
        if gciAbs(1)-frameStart(end) > winsz %20e-3*fs
            %frameStart(end+1)=segStart;
            %frameStart(end+1)=round((gciAbs(1)+frameStart(end))/2);
            newFrSt=(frameStart(end)+winsz):winsz:(gciAbs(1));
            if gciAbs(1)-newFrSt(end)<2.5e-3*fs
                newFrSt(end)=[];
            end
            frameStart=[frameStart newFrSt];
        end
            
        frameStart=[frameStart gciAbs]; %         
    end
    
    % Next unvoiced frame
    if frameStart(end) < Nlen-2*winsz
        segStart=frameStart(end)+winsz;
        if ii<length(iT)  % Check for end of speech
            segEnd=iT(ii+1);
        else
            segEnd=Nlen;
        end
        
        
        frameStart=[frameStart segStart:winsz:segEnd];
    end
end

figure(22); clf;
plot(1000*diff(frameStart)/fs); hold on;
xlabel('Samples')
ylabel('Time [ms]')
title('Length of windows (diff(frameStart))')

figure(23); clf;
plot(numSmall);



    
    
