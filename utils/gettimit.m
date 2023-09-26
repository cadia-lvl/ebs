function [s,fs,wrd,phn,txt]=gettimit(f,m)
%  Inputs:   F       file name ['TEST/DR1/FAKS0/SA1.WAV']
%            M       mode string:
%                    'n'  Normalize speech level to 0 dB using ITU P.56
%                    'v'  Trim leading and ending silences
%
% Outputs:   S(:,1)  speech waveform
%           FS       sample frequency (always 16000)
%          WRD(nw,2} each row is of the form {[t1 t2] 'word'} where t1 and t2 give the
%                    start end end times in seconds and 'wor' gives the text of the word
%          PHN{np,2} each row is of the form {[t1 t2] 'phon'} where t1 and t2 give the
%                    start end end times in seconds and 'wor' gives the text of the phoneme
%          TXT(1,:)  Text of the sentence
%
% TIMIT contains 630 talkers each saying 10 sentences: 2 SA, 5 SX and 2 SI
% sentences. The sentence types are: SA = dialect identification, SX =
% phonetically compact to provide a good coverage of phone pairs, SI =
% phonetically diverse from existing text sources. All speakers read the
% same two SA sentences (normally omitted from evaluations), seven speakers
% read each SX sentence and each SI sentence is unique.
% The traning set includes 462 speakers. The core test set contains 24 speakers
% (2 male and 1 female speaker from each of 8 dialect regions) while the full
% test set includes 168 speakers. No SI or SX sentences are common to the test
% and training sets.
%
% see also:  timitfiles
persistent dd

if isempty(dd)
    par=projParam();
    dd=par.pth.speechpth;
end
if nargin==0 || ~numel(f)
    f='TEST/DR1/FAKS0/SA1.WAV'; % default file
end
f(f=='\')='/';  % convert the delimiters
if f(end)=='/'
    f(end)=[];
end
nd=sum(f=='/'); % number of delimiters
if nd==2
    f=[f '/SA1.WAV'];
end
if nargin>1 && length(m)>0
    mv=any(m=='v'); % apply VAD to trim silences
    mn=any(m=='n'); % normalize speech level to 0 dB
else
    mv=false;
    mn=false;
end
nag=max(nargout,4*mv);
switch nag
    case {1,2}
        [s,fs]=v_readsph([dd f],'wt');
    case {0,3}
        [s,fs,wrd]=v_readsph([dd f],'wt');
    case 4
        [s,fs,wrd,phn]=v_readsph([dd f],'wt');
    case 5
        [s,fs,wrd,phn]=v_readsph([dd f],'wt');
        nw=size(wrd,1);
        if ~nw
            txt='';
        else
            nc=0;
            for i=1:nw
                nc=nc+numel(wrd{i,2});
            end
            txt=repmat(' ',1,nc+nw-1);
            j=1;
            for i=1:nw
                ti=wrd{i,2};
                k=j+numel(ti)+1;
                txt(j:k-2)=ti;
                j=k;
            end
        end
end
if mv
    if strcmp(phn{end,2},'h#')                      % if the final phone is silence
        s(round(phn{end,1}(1)*fs+1):end)=[];       % trim the speech
        phn(end,:)=[];                              % and remove the silence phone
    end
    if strcmp(phn{1,2},'h#')                        % if the first phone is silence
        ntrim=round(phn{1,1}(2)*fs);                % samples to trim
        s(1:ntrim)=[];                             % trim the speech
        phn(1,:)=[];                                % and remove the silence phone
        jmax=size(phn,1);                           % number of remaining phones
        for j=1:jmax                                % loop through the phones
            phn{j,1}=(round(phn{j,1}*fs)-ntrim)/fs; % update the phoneme positions
        end
        jmax=size(wrd,1);                           % number of words
        for j=1:jmax                                % loop through the words
            wrd{j,1}=(round(wrd{j,1}*fs)-ntrim)/fs; % update the word positions
        end
    end
end
if nargin>1 && ~isempty(m)
    if any(m=='n')                                  % normalize speech level to 0 dB
        s=v_activlev(s,fs,'n');
    end
end
if ~nargout
    % Plot a spectrogram with the word transcription
    spgrambw(s,fs,'Jwcpta',[],[],[],[],wrd);
    title(['TIMIT/' f]);
end