% Print list of TIMIT phones
dofreq=1;   % 0=omit occurrence frequencies, 1=find occurrence frequencies
timit=gettimitpath;                         % get path to timit subfolder of timit database
nphn=67; % number of TIMIT phones (including four that are not in transcriptions: 'cl','vcl','1','2')
[st65,vt65]=w_phoncode('t',1:nphn); % get TIMIT65 strings
u65=w_phoncode('tU',1:nphn); % get TIMIT65 strings
[st48,vt48]=w_phoncode('tf',st65); % get TIMIT48 strings
[st39,vt39]=w_phoncode('tF',st65); % get TIMIT48 strings
cnts=zeros(nphn,2); %iniialize counts
if dofreq
    tfarg={'n' 't'}; % argument for timitfiles call
    for it=1:2 % loop for train and test
        tf=timitfiles(tfarg{it},0,timit); % get list of training files
        ntf=length(tf);
        for i=1:ntf
            v_finishat([i 1 ntf]);
            fnam=[timit tf{i}]; % get the .wav filename
            fnam(end-2:end)='PHN'; % convert to the phonetic transcription filename
            fid=fopen(fnam,'r');
            if fid<0
                error('Cannot open %s',[timit tf{i}]);
            end
            cphn = textscan(fid,'%d%d%s');
            fclose(fid);
            phnlist=string(cphn{3}); % convert phone list to strings
            [sph,vph]=w_phoncode('t',phnlist); % get TIMIT65 strings
            cnts(:,it)=cnts(:,it)+sparse(vph(:,1),1,1,nphn,1); % add into phone counts
        end
    end
end
%
fprintf(' N TIMIT TIM48 TIM39 Class Place Features  Train  Test Unicode\n');
for i=1:nphn
    if strcmp(st48(i),st65(i))
        st48(i)="";
    end
    if strcmp(st39(i),st65(i))
        st39(i)="";
    end
    fprintf('%2d%5s%6s%6s%5d%6d   %08x%7d%6d%5s\n',i,st65(i),st48(i),st39(i),vt65(1,i,3:5),cnts(i,:),u65(i));
end