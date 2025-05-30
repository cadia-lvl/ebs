function [sp,fs,timeMarkers] = getSpeechAndTimeMarkers(spfile,fs,fmt)

% spfile -  Name of speech file to be processed
% tmpDir -  Directory (absolute path) to store temporary wav files 

%praatScriptDir='/Users/jg/Work/PraatTools';
praatScriptDir='/Users/jg/GoogleDrive/Work/Software/Praat/';
%fs=20000;
if strcmp(fmt,'T')
    [spold,fsold]=readsph(spfile);
    if fs==0
        fs=fsold;
    end
    sp=resample(spold,fs,fsold);
    phnfile=[spfile(1:end-3) 'PHN'];
    phn=extractPhn(phnfile);
    dvu=diff(phn.vus==3);
    timeMarkers(:,1)=find(dvu==1)/fs;
    timeMarkers(:,2)=find(dvu==-1)/fs;
else
    [spold,fsold] = readwav(spfile);
    
    if fs==0
        fs=fsold;
    end
    
    tmpInFile=[praatScriptDir '/tmp0.wav'];
    try
        sp=resample(spold(:,2),fs,fsold);
        writewav(sp,fs,tmpInFile);
        timeMarkers = vuvExtract(tmpInFile,praatScriptDir);
        delete(tmpInFile);
    catch
        sp=resample(spold(:,1),fs,fsold);
        writewav(sp,fs,tmpInFile);
        timeMarkers = vuvExtract(tmpInFile,praatScriptDir);
        delete(tmpInFile);
    end
end