function [infoStruct, gci] = fftAnalysisOfLongSpeechSignal(sp, fs, timeMarkers)

iT=round(timeMarkers(:,1)*fs);
eT=round(timeMarkers(:,2)*fs);

%Extend (dilate) the voiced segment by 30 ms either side:
nExt= round(30e-3*fs);

gci=[];
goi=[];
%vsf=[];
%segIx=[];
for ii=100%1:length(iT)
    % Make sure the segments to index out of the signal
    segStart = max(iT(ii)-nExt,1);
    segEnd   =  min(eT(ii)+nExt,length(sp));
    
    if (segEnd-segStart)< fs*50e-3   %Ignore small (or even negative) segments
        continue;
    end
    
    ignoreGOIfromYaga=1;
    if ignoreGOIfromYaga %
        gciSeg = dypsagoi(sp(segStart:segEnd),fs);
        gci(end+1:end+length(gciSeg)) = gciSeg+segStart-1;
    else
        
        [gciSeg, goiSeg] = dypsagoi(sp(segStart:segEnd),fs);
        %gci = [gci, gciSeg+segStart-1];
        %goi = [goi, goiSeg+segStart-1];
        
        gci(end+1:end+length(gciSeg)) = gciSeg+segStart-1;
        goi(end+1:end+length(goiSeg)) = goiSeg+segStart-1;
        
        infoStruct(ii).goiSeg = goiSeg;
    end;
    [spseg] = pitchSynchronousFreqAnalysis(sp(segStart:segEnd), fs, gciSeg);


    infoStruct(ii).gciSeg = gciSeg;
    infoStruct(ii).segStart = segStart;
    infoStruct(ii).segEnd = segEnd;
    
end;





    
    
