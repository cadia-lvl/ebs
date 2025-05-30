function [vsfeatures, gci, u, udash] = vsProcessLongSpeechSignal(sp, fs, timeMarkers)



iT=round(timeMarkers(:,1)*fs);
eT=round(timeMarkers(:,2)*fs);

[udash, ~, ~, u] = iaif(sp,fs);

%Extend (dilate) the voiced segment by 30 ms either side:
nExt= round(30e-3*fs);

gci=[];
goi=[];
%vsf=[];
%segIx=[];
for ii=1:length(iT)
    %Make sure the segments to index out of the signal
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
        
        vsfeatures(ii).goiSeg = goiSeg;
    end;
    
    [mfdr, cq, pa, naq, f0, h1h2, hrf] = extractVoiceFeatures(u(segStart:segEnd), fs, gciSeg);
    
    
    vsfeatures(ii).mfdr =mfdr;
    vsfeatures(ii).cq =cq;
    vsfeatures(ii).pa = pa;
    vsfeatures(ii).naq = naq;
    vsfeatures(ii).f0 = f0;
    vsfeatures(ii).h1h2 =h1h2;
    vsfeatures(ii).hrf = hrf;
    vsfeatures(ii).gciSeg = gciSeg;
    vsfeatures(ii).segStart = segStart;
    vsfeatures(ii).segEnd = segEnd;
    
    %vsf(:,end+1:end+length(mfdr),:) = [mfdr; cq; pa; naq; f0; h1h2; hrf; gciSeg];
    %segIx(end+1:end+length(mfdr))=ii;
    
end;





    
    
