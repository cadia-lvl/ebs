function [sout,vout]=w_phoncode(mode,sin,msk)
% Usage: (1) sout=w_phoncode('tU',"tcl ch ow kcl k");    % convert "choke" from TIMIT to Unicode
%
%        (2) [st,vt]=w_phoncode('t',1:67);                      % get all TIMIT strings
%            u=w_phoncode('tU',1:67);                           % convert to Unicode
%            fprintf(' N TIMIT Class Place Features  Uni\n');   % print list of TIMIT phones
%            for i=1:67                                         % loop through all 67 phones
%              fprintf('%2d%5s%5d%6d   %08x%5s\n',i,st(i),vt(1,i,3:5),u(i));
%            end
%
%  Inputs: mode     Character array giving the following options:
%                       't'     format is TIMIT: 65 phones + 2 stress marks
%                       's'     format is SILIPA93: 255 phones
%                       'u'     sin input gives IPA numbers or Unicode phone strings instead of the numbers/symbols from the chosen format
%                       'U'     sout output should be Unicode characters instead of those from the chosen format
%                       'f'     fold TIMIT allophones to reduce to 48 phones + 2 stress marks
%                                   phone-numbers 8,10,12,15,25,40,42,52-59,63,65 are unused
%                       'F'     fold TIMIT allophones to reduce TIMIT to 39 phones + 2 stress marks
%                                   phone-numbers  8,10,12,13,15,23,25,29,37,38,40,42,44,52-63,65 are unused
%                       'j'     join sout phones into a single string with space delimiters
%                       'J'     join sout phones into a single string without delimiters
%          sin      Input phone sequence as string or character array (trailing blanks are removed) or as vector of phone numbers
%          msk      Selects a subsets of Distinctive Features to be included in the vout(:,5) bit mask
%
% Outputs: sout     Output phone sequence as string or character array
%                       As a special case, if w_phoncode is called with no input arguments
%                       then sout(30,3) lists the codes shows below for class, place and features.
%          vout     Numerical array with phone info: [phone-number IPAnum Class Place Features]
%                       Notes:  (a) vout(:,1) equals the phone number of the chosen format
%                               (b) vout(:,2) equals IPA2*1000+IPA1 where IPA1 amd IPA2 are the IPA number of the phone or phone pair
%
%                       vout(:,3) Class values:
%                               1=Silence,      2=Nonspeech,    3=Plosive,              4=Nasal,                    5=Trill,
%                               6=Flap,         7=Percussive,   8=Sibilant-Fricative,   9=Non-sibilant Fricative,   10=Lateral-fricative,
%                               11=Affricate,   12=Approximant, 13=Lateral-approximant, 14=Vowel,                   15=Dipthong,
%                               16=Implosive,   17=Click,       18=Ejective,            19=Closure,                 20=Voicing,
%                               21=Diacritic,   22=Stress,      23=Tone,                24=Prosody,                 25=Suprasegmental,
%                               26=Delimitation
%
%                       vout(:,4) Place values:
%                               1=Bilabial,     2=Bidental,         3=Labiodental,      4=Interdental,      5=Dental,
%                               6=Alveolar,     7=Palatoalveolar,   8=Retroflex,        9=Alveolopalatal,   10=Palatal,
%                               11=Labiovelar,  12=Velar,           13=Velopharangeal,  14=Uvular,          15=Pharyngeal,
%                               16=Epiglottal,  17=Glottal
%
%                       vout(:,5) Distinctive Feature default bit positions (use the mask input to select only a subset):
%                               0=silence,          1=syllabic,         2=consonantal,          3=sonorant,             4=continuant,
%                               5=delayed release,  6=strident,         7=nasal,                8=lateral,              9=tap-or-flap,
%                               10=trill,           11=voice,           12=spread glottis,      13=constricted glottis, 14=implosive,
%                               15=pharangeal,      16=dorsal,          17=distributed,         18=anterior,            19=coronal,
%                               20=labiobidental,   21=labial,          22=round,               23=low,                 24=mid,
%                               25=high,            26=back,            27=front,               28=tense,               29=dipthong
%
% Notes:
%       (1) no TIMIT transcriptions include the stress marks or the two phones 60="cl" and 61="vcl"
%
% Bugs/suggestions:
%       (1) Other formats: 'm'=SAMPA, 'x'=X-SAMPA, 'a'=ARPABET, 'A'=case-sensitive-ARPABET,
%                          'S'=SILIPA90, 'k'=IPA-SIL-keyboard, 'p'=PRAAT, 'K'=Kirshenbaum-usenet, 'l'=LJSpeech
%       (2) Optionally simplify an input phone string
%       (3) Default should not be TIMIT but Unicode
%
persistent phonv phonstr ipad fmts fmtc unid key
if isempty(phonstr)
    key = [ ...
        "Silence"                   "Bilabial"          "silence"               ;% Word  1, Bit    0
        "Nonspeech"                 "Bidental"          "syllabic"              ;% Word  2, Bit    1
        "Plosive"                   "Labiodental"       "consonantal"           ;% Word  3, Bit    2
        "Nasal"                     "Interdental"       "sonorant"              ;% Word  4, Bit    3
        "Trill"                     "Dental"            "continuant"            ;% Word  5, Bit    4
        "Flap"                      "Alveolar"          "delayed release"       ;% Word  6, Bit    5
        "Percussive"                "Palatoalveolar"    "strident"              ;% Word  7, Bit    6
        "Sibilant-Fricative"        "Retroflex"         "nasal"                 ;% Word  8, Bit    7
        "Non-sibilant Fricative"    "Alveolopalatal"    "lateral"               ;% Word  9, Bit    8
        "Lateral-fricative"         "Palatal"           "tap-or-flap"           ;% Word 10, Bit    9
        "Affricate"                 "Labiovelar"        "trill"                 ;% Word 11, Bit   10
        "Approximant"               "Velar"             "voice"                 ;% Word 12, Bit   11
        "Lateral-approximant"       "Velopharangeal"    "spread glottis"        ;% Word 13, Bit   12
        "Vowel"                     "Uvular"            "constricted glottis"   ;% Word 14, Bit   13
        "Dipthong"                  "Pharyngeal"        "implosive"             ;% Word 15, Bit   14
        "Implosive"                 "Epiglottal"        "pharangeal"            ;% Word 16, Bit   15
        "Click"                     "Glottal"           "dorsal"                ;% Word 17, Bit   16
        "Ejective"                  ""                  "distributed"           ;% Word 18, Bit   17
        "Closure"                   ""                  "anterior"              ;% Word 19, Bit   18
        "Voicing"                   ""                  "coronal"               ;% Word 20, Bit   19
        "Diacritic"                 ""                  "labiobidental"         ;% Word 21, Bit   20
        "Stress"                    ""                  "labial"                ;% Word 22, Bit   21
        "Tone"                      ""                  "round"                 ;% Word 23, Bit   22
        "Prosody"                   ""                  "low"                   ;% Word 24, Bit   23
        "Suprasegmental"            ""                  "mid"                   ;% Word 25, Bit   24
        "Delimitation"              ""                  "high"                  ;% Word 26, Bit   25
        ""                          ""                  "back"                  ;% Word 27, Bit   26
        ""                          ""                  "front"                 ;% Word 28, Bit   27
        ""                          ""                  "tense"                 ;% Word 29, Bit   28
        ""                          ""                  "dipthong"   ]          ;% Word 30, Bit   29
    load([mfilename('fullpath') '.mat'],'ps');   % load .mat file containing data
    fmtc=['tl'; 'sx']; % one row per format: [id-char case-insensitive]
    nfmt=size(fmtc,1); % number of formats
    fmts=cell(nfmt,5);
    fmts(:,1)={'TIMIT';'SILIPA93'};
    %     each phone has multiple numbers: (a) phone num = index into list, (b) phone id = index in original database,
    %                                      (c) [optional] IPA num = IPA1 + 1000*IPA2 where IPA1 and IPA2 are one or two IPA numbes
    %                                      (d) [optional;] one or more format numbers, e.g. TIMIT number
    nphon=length(ps.Unicode1);              % number of phones in list
    seqlist=(1:nphon)';                     % list of phone numbers
    phonstr=deblank(string(char([ps.Unicode1 ps.Unicode2 ps.Unicode3]))); % convert unicode to string array
    % sort out PhoneID mappings
    phid=ps.PhoneID;
    maxid=max(phid);                    % find maximum PhoneID
    phseq=zeros(maxid,1);               % space for PhoneID->phseq
    phseq(phid)=1:length(phid);         % mapping from PhoneID to sequence number (0 if none)
    % sort out preferred alternatives
    phlist=ps.NPreferred;               % list of preferred phones
    phlist(phlist>maxid)=NaN;           % remove any PhonIDs > maxid
    phpref=zeros(nphon,1);                 % initialize each to zero
    phpref(~isnan(phlist))=phseq(phlist(~isnan(phlist))); % replace PhonIDs by sequence indexes (0 if no preferred alternative)
    % sort out equivalents
    phlist=ps.NEquiv;                   % list of equivalent phones
    phlist(phlist>maxid)=NaN;           % remove any PhonIDs > maxid
    phequiv=zeros(nphon,1);             % initialize each to point to itself
    phequiv(~isnan(phlist))=phseq(phlist(~isnan(phlist))); % replace PhonIDs by sequence indexes (0 if no preferred alternative)
    % create numerical phon array
    phonv=[ps.IPAn ps.NClass ps.NPlace ps.Feats phpref phequiv]; % table of phone information indexed by phone number
    phonv(isnan(phonv))=0; % replace NaNs by 0s
    ipaseq=sortrows([phonv(:,[5 1]) (1:nphon)'],'descend'); % put ipa codes with preferences first so they are overwritten
    ipad=dictionary(ipaseq(:,2),ipaseq(:,3)); % create dictionary: IPA number to phone number
    unid=dictionary(phonstr,seqlist); % create dictionary: Unicode string to phone number
    % sort out format mappings
    for i=1:nfmt
        fmti=fmts{i,1}; % format name
        fmtin=[fmti 'n']; % format name + 'n' lists format numbers
        psfmti=deblank(string(ps.(fmti)));
        psfmtin=ps.([fmti 'n']); % list of format numbers
        mk=strlength(psfmti)>0;                    % mask for defined phones
        timn=psfmtin(mk);               % primary format numbers (not necessarily in numerical order)
        timmax=max(timn);                       % highest format code
        if fmtc(i,1)=='t'
            timitv=zeros(timmax,3); % create space for mapping from format num to phone num and to reduced sets of phone num
            timitv(timn,:)=[seqlist(mk) ps.TIMIT48n(mk) ps.TIMIT39n(mk)];
        else
            timitv=zeros(timmax,1); % create space for mapping from format num to phone num
            timitv(timn,:)=seqlist(mk); % create mapping vector
        end
        timits=repmat("",timmax,1);
        timits(timn)=psfmti(mk);                   % valid format strings
        fmts{i,2}=timits;
        fmts{i,3}=timitv;
        fmts{i,4}=dictionary(timits(timn),timn); % map format string to format num
        mk=~isnan(psfmtin); % all valid format numbers
        fmts{i,5}=dictionary(seqlist(mk),psfmtin(mk)); % map phone num to format num
    end
end
if nargin==0
    sout=key;
    vout=[];
else
    ifmt=1; % TIMIT is default format
    if(any(mode=='t')) % this could be cleaned up, e.g. find(any(repmat(mode',1,size(fmtc,1))==repmat(fmtc,length(mode),1),1))
        ifmt=1;
    elseif(any(mode=='s'))
        ifmt=2;
    % elseif(any(mode=='u') && any(mode=='U')) % should cope better with pure unicode
    %     ifmt=0; % purely unicode
    end
    nphfmt=size(fmts{ifmt,3},1); % number of phones in this format
    if nargin<2 || isempty(sin)             % there is no sin input, so list all phones of requested type
        nsin=size(fmts{ifmt,3},1);      % number of timit codes
        seqs=fmts{ifmt,3}(:,1);         % and corresponding phone numbers
    else                                    % sin input given
        if ischar(sin) || iscell(sin)
            sin=string(sin);                % convert char input to string
        end
        if isstring(sin)                    % if sin is char or string
            if numel(sin)==1
                sins=split(deblank(sin));   % if only one element, split it into phones
            else
                sins=deblank(sin(:));       % else make it a column vector
            end
            nsin=size(sins,1);              % number of input phones
            seqs=zeros(nsin,1); % space for phone numbers
            if any(mode=='u')               % sin contains unicode strings
                mk=isKey(unid,sins); % mask for strings in dictionary
                seqs(mk)=unid(sins(mk)); % convert unicode strings to phone numbers
            else
                if fmtc(ifmt,2)=='l'
                    sins=lower(sins);           % TIMIT is case-insensitive
                end
                mk=isKey(fmts{ifmt,4},sins); % mask for strings in dictionary
                timix=zeros(nsin,1);
                timix(mk)=fmts{ifmt,4}(sins(mk)); % convert strings to format numbers
                if fmtc(ifmt,1)=='t' % for TIMIT only
                    if any(mode=='F')
                        timix(mk)=fmts{ifmt,3}(timix(mk),3); % reduce to 39 phone set
                    elseif any(mode=='f')
                        timix(mk)=fmts{ifmt,3}(timix(mk),2); % reduce to 48 phone set
                    end
                end
                seqs(mk)=fmts{ifmt,3}(timix(mk),1); % convert format numbers to phone numbers
            end
        else                                % sin is numeric
            sins=sin(:);                    % make a column vector
            nsin=size(sins,1);              % number of input phones
            seqs=zeros(nsin,1);
            if any(mode=='u')               % sin contains IPA numbers
                mk=isKey(ipad,sins);
                seqs(mk)=ipad(sins(mk));    % convert to phone numbers
            else % sin contains format numbers
                mk=(sins>0) & (sins<=nphfmt); % mask for valid codes
                seqs(mk)=fmts{ifmt,3}(sins(mk),1);   % convert to phone numbers
            end
        end
    end
    vout=zeros(nsin,5);
    mk=seqs>0;
if(any(mode=='u') && any(mode=='U'))
    vout(mk,1)=seqs(mk); % use internal IPA phone numbers
else
    vout(mk,1)=fmts{ifmt,5}(seqs(mk)); % convert phone numbers to format numbers for column 1
end
    vout(mk,2:5)=phonv(seqs(mk),1:4);
    if nargin>2 && ~isempty(msk)
        vout(:,5)=rem(floor(vout(:,5)*pow2(-msk(:)')),2)*pow2(0:length(msk)-1)'; % pick out the requested bits
    end
    sout=repmat("",nsin,1);
    if any(mode=='U')
        sout(mk)=phonstr(seqs(mk));
    else
        mk=vout(:,1)>0;
        sout(mk)=fmts{ifmt,2}(vout(mk,1));
    end
    % reshape vout and sout if necessary
    if nargin>1
        sz=size(sin);
        if (sz(2)>1 || length(sz)>2)
            sout=reshape(sout,size(sin));
            vout=reshape(vout,[size(sin) 5]);
        end
    end
    if any(lower(mode)=='j')
        if any(mode=='j')
            sout=join(sout(:)); % 'j': reassemble as a single input string with space delimiter
        else
            sout=join(sout(:),''); % 'J': reassemble as a single input string without a delimiter
        end
    end
end

