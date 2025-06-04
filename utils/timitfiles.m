function [tf,ty,tk,s,fs,wrd,phn,txt]=timitfiles(m,n,timit)
%  Inputs:  m  string specifying the type of TIMIT file to extract:
%                'c'  restrict to core test set
%                'A'  include SA sentences (normally excluded)
%                'a'  restrict to SA sentences (dialect ID, each read by all 630 speakers)
%                'i'  restrict to SI sentences (diverse, each read by 1 speaker)
%                'x'  restrict to SX sentences (phonetically compact, each read by 7 speakers)
%                'n'  restrict to train
%                't'  restrict to test
%                'm'  restrict to male
%                'f'  restrict fo female
%                '1'-'8' restrict to dialect region
%                'z'  don't use cell arrays for outputs tf, s, wrd, phn and txt when only one file is output
%                'v'  remove silent segemnts at the beginning and end
%                'p'  normalize level with ITU P.56
%                'j'  Join all utterances together in s, wrd, phn and txt outputs
%              alternatively, m can be a case-insensitive regular expression specifying the wanted
%              files (which must include '/'). Default is 't(est|rain)/dr./...../s(i|x).*\.wav$'
%           n  choose n files randomly from matches (or all if n=0)
%       timit  base folder for the TIMIT files with or without a final '/' (parent of TEST and TRAIN)
%
% Outputs: tf{n,1}   cell vector containing file names e.g. 'TEST/DR1/FAKS0/SI1573.WAV'
%          ty(n,12)  char array with each row containing 'trgxxxnsyyyy' where
%                       t = T or N  for Test and Train sets
%                       r = 1 to 8 for dialect region
%                       g = F or M for gender of talker
%                    xxxn = speaker initials ending with disambiguating digit (usually 0)
%                       s = A, I or X for type of sentence
%                    yyyy = token number
%          tk(:,5)   sorted char array with a list of unique talkers in the form 'gxxxn' (see above)
%          s{n,1}    speech waveforms (as cell vector of column vectors)
%          fs        sample frequency (always 16000)
%          wrd{n,1}  Each entry is a cell array, {nw,2}, where each row is of the form {[t1 t2] 'word'} where t1 and t2 give the
%                    start end end times in seconds and 'word' gives the text of the word
%          phn{n,1}  Each entry is a cell array, {np,2}, where each row is of the form {[t1 t2] 'phon'} where t1 and t2 give the
%                    start end end times in seconds and 'wor' gives the text of the phoneme
%          txt{n,1}  Text of the sentences as cell vector
%
% The first call to timitfiles() willl be slow because the routine reads
% the entire TIMIT directory but subsequent calls will be much quicker.
persistent dd corem coref flregs ddt tff tyy
if isempty(corem)
    corem='DAB0|WBT0|TAS1|WEW0|JMP0|LNT0|LLL0|TLS0|BPM0|KLT0|CMJ0|JDH0|GRT0|NJM0|JLN0|PAM0';
    coref='ELC0|PAS0|PKT0|JLM0|NLP0|MGD0|DHC0|MLD0';
    flregs={'(i|x)' 'x' 'i' '(i|x)' 'a' '(a|x)' '(a|i)' '(a|i|x)'};
    cnam=computername;
    if strcmp(cnam,'ee-dmblap2')
        dd='C:/OneDriveImperial/OneDrive - Imperial College London/work/data/speech/timit/timit/';
    else
        dd='D:/OneDrive - Imperial College London/work/data/speech/timit/timit/';
    end
    ddt='';
end
if nargin<3 || isempty(timit)
    timit=dd;
end
if exist(fullfile(timit,'test'),'dir')+exist(fullfile(timit,'train'),'dir')==0
    if exist(fullfile(timit,'timit','test'),'dir')+exist(fullfile(timit,'timit','train'),'dir')==0
        error('Cannot find subfolders TEST or TRAIN within TIMIT folder:\n    %s',timit);
    else
        error('Cannot find subfolders TEST or TRAIN within TIMIT folder:\n    %s\nPlease set TIMIT folder instead to:\n    %s',timit,fullfile(timit,'timit'));
    end
end
if ~strcmpi(timit,ddt)
    tff=v_regexfiles('.wav$',timit,'r');
    if isempty(tff)
        error('No speech files found at %s',timit);
    end
    ddt=timit;
    %     tyc=cellfun(@(x) regexp(x,'^[^/]*(.)/[^/]*(.)/(.)[^/]*/.(.)','tokens'),tff); % extract group, region, gender, type
    tyc=cellfun(@(x) regexp(x,'^[^/]*(.)/[^/]*(.)/([^/]*)/.(.)([^.]*)','tokens'),tff); % extract group, region, gender, type
    tyy=repmat('?',length(tyc),12); % reserve space
    for i=1:length(tyc)
        tyy(i,:)=sprintf('%s%04d',cell2mat(tyc{i}(1:4)),str2num(tyc{i}{5}));
    end
end
if nargin<1
    m='';
end
if isempty(m)
    regex='t(est|rain)/dr./...../s(i|x).*\.wav$';
elseif any(m=='/')
    regex=m;
else
    % select test/train subset
    if any(m=='c') || any(m=='t') && ~any(m=='n')
        sbreg='test';
    elseif any(m=='n') && ~any(m=='t')
        sbreg='train';
    else
        sbreg='t(est|rain)';
    end
    % select dialect region
    drmode=m>='1' & m<='8';
    if any(drmode)
        drreg=m(drmode);
        drreg=['dr(' reshape([drreg; repmat('|',1,length(drreg)-1) ')'],1,[])];
    else
        drreg='dr.';
    end
    % select speaker
    if any(m=='c')
        if any(m=='m') && ~any(m=='f')
            spreg=['m(' corem ')'];
        elseif any(m=='f') && ~any(m=='m')
            spreg=['f(' coref ')'];
        else
            spreg=['(m|f)(' corem '|' coref ')'];
        end
    elseif any(m=='m') && ~any(m=='f')
        spreg=['m....'];
    elseif any(m=='f') && ~any(m=='m')
        spreg=['f....'];
    else
        spreg='.....';
    end
    spmode=m=='m' | m=='f';
    % select file
    flreg=4*any(m=='a')+2*any(m=='i')+any(m=='x')+1;
    if flreg==1 && any(m=='A')
        flreg='(a|i|x)';
    else
        flreg=flregs{flreg};
    end
    regex=[sbreg '/' drreg '/' spreg '/s' flreg '.*\.wav$'];
end
mk=~cellfun('isempty',regexpi(tff,regex));
tf=tff(mk);
ty=tyy(mk,:);
ntf=size(ty,1);
if nargin>1 && n>0 && n<ntf
    ix=v_randiscr(ones(ntf,1),-n); % pick n files at random
    tf=tf(ix);
    ty=ty(ix,:);
    ntf=size(ty,1);
end
if nargout>2 || ~nargout
    tk=unique(ty(:,3:7),'rows'); % get list of talkers
    if nargout>3 || ~nargout
        s=cell(ntf,1);
        wrd=cell(ntf,1);
        phn=cell(ntf,1);
        txt=cell(ntf,1);
        if any(m=='p')
            gtm='n';        % normalize option for gettimit
        else
            gtm='';
        end
        if any(m=='v')      % perform VAD and remove top & tail silences
            gtm=[gtm 'v'];
        end
        for i=1:ntf
            [s{i},fs,wrd{i},phn{i},txt{i}]=gettimit(tf{i},gtm);
        end
    end
end
if ntf>1 && any(m=='j') % join files together
    ns=size(s{1},1); % number of speech samples so far
    nwrd=size(wrd{1},1); % number of words so far
    nphn=size(phn{1},1); % number of phones so far
    for i=2:ntf
        s{1}=cat(1,s{1},s{i}); % append speech waveform
        wrd{1}=cat(1,wrd{1},wrd{i}); % append speech waveform
        phn{1}=cat(1,phn{1},phn{i}); % append speech waveform
        nmx=size(wrd{1},1); % new number of words
        for j=nwrd+1:nmx
            wrd{1}{j,1}=(round(wrd{1}{j,1}*fs)+ns)/fs;
        end
        nwrd=nmx; % update number of words
        nmx=size(phn{1},1); % new number of phones
        for j=nphn+1:nmx
            phn{1}{j,1}= (round(phn{1}{j,1}*fs)+ns)/fs;
        end
        nphn=nmx; % update number of words
        ns=size(s{1},1); % update number of speech samples
        txt{1}=[txt{1} ' ' txt{i}];
    end
    s(2:ntf)=[]; % remov obsolete cells
    wrd(2:ntf)=[];
    phn(2:ntf)=[];
    txt(2:ntf)=[];
    ntf=1;
end
if ~nargout
    % Plot a spectrogram with the word transcription
    spgrambw(s{1},fs,'Jwcpta',[],[],[],[],wrd{1});
elseif ntf==1 && any(m=='z') % omit cell arrays if only one output
    tf=tf{1};
    s=s{1};
    wrd=wrd{1};
    phn=phn{1};
    txt=txt{1};
end

