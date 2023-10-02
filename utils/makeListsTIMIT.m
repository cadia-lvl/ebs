function makeListsTIMIT(par,grepstr)

%   Makes a file list for a specific database specified in par
%
%   Inputs:
%       par                     A structure containing with information
%                               about paths such as output files etc,
%               .speechpth      Path to top of the TIMIT tree (the top
%                               TIMIT directory)
%               .frm            Format ending of the wave files
%               .fetpth         Path to the top of the feature directory
%               .splistname     The name of the speech file list, e.g. splist.scp
%               .fetlistname    The name of the feature file list, e.g.fetlist.scp
%               .phlistname     The name of the phone file list, e.g.phlist.scp
%               .listopt         'sfp';     %(s = make speech file list, f = make file output list, f = make phn list)
%               .fetext         Extention of the feature file, e.g. 'mat';
%       grepstr          (optional) grep the list for this string, useful
%                        for specifying a subset of the data
%
%   References:
%
if nargin<2
    grepstr='';
end
%Short-hands
listopt=par.listopt;
frm=par.frm;
splistname=par.splistname;
fetlistname=par.fetlistname;
phlistname=par.phlistname;
fetext=par.fetext;
fetpth=par.fetpth;
if ~numel(fetpth)
    fetpth='.';                         % set fetpth to the current directory
end
if ispc
    tok=v_regexfiles([grepstr '.*\.' frm '$'],par.speechpth,'r');   % recursively find files matching grepstr
    nfrm=length(frm);
    if any(listopt=='s')                                            % option 's': Create a list of speech files
        fid=fopen(splistname,'wt');
        for i=1:length(tok)
            fprintf(fid,'%s\n',fullfile(par.speechpth,tok{i}));
        end
        fclose(fid);
    end
    if any(listopt=='p')                                            % option 'p': Create a list of phoneme files
        fid=fopen(phlistname,'wt');
        for i=1:length(tok)
            toki=tok{i};
            fprintf(fid,'%s\n',fullfile(par.speechpth,[toki(1:end-nfrm) 'PHN']));
        end
        fclose(fid);
    end;
    if any(listopt=='f')                                            % option 'f': Create a list of feature files and create directory structure
        maked=1;                                                    % default is to make necessary directory structure
        if ~exist(fetpth,'dir') || length(dir(fetpth))>2
            disp(['Not creating path structure since ' fetpth ' either doesn''t exist or is not empty']);
            maked=0;                                                % do not make directory structure
        end;
        fid=fopen(fetlistname,'wt');
        prevpath='///';                                             % initialize to impossible path
        for i=1:length(tok)
            toki=tok{i};
            if maked
                path=fileparts(toki);                               % get path of speech file
                if ~strcmpi(path,prevpath)                          % skip directory creation if the same as the previous file
                    fullpath=fullfile(fetpth,path);                 % preappend the feature path
                    if ~exist(fullpath,'dir');
                        mkdir(fullpath);                            % create the directory
                    end
                end
                prevpath=path;                                      % save this path for next time around
            end
            fprintf(fid,'%s\n',fullfile(fetpth,[toki(1:end-nfrm) fetext]));
        end
        fclose(fid);
    end;
else
    if any(listopt=='s')                                                    % Make splist
        eval(['!(cd ' par.speechpth '; find * -type f -name "*.' frm '") > tmp1.scp;']);    % find all speech files
        if ~isempty(grepstr)
            eval(['!egrep "' grepstr '" tmp1.scp > tmp2.scp']);             % only keep those that match grepstr
            eval(['!rm -f tmp1.scp']);
        else
            eval(['!mv tmp1.scp tmp2.scp']);                                % keep all files if grepstr is empty
        end;
        eval(['!sed "s:^:' par.speechpth '/:" tmp2.scp > ' splistname]);    % prefix with the path to the top directory
        eval(['!rm -f tmp2.scp']);                                          % delete the temporary file
    end;

    if any(listopt=='p')
        eval(['!sed "s:' frm ':PHN:" ' splistname ' > ' phlistname]);       % replace the file extension with PHN
    end;

    if any(listopt=='f')
        dd=dir(fetpth);
        if length(dd)>2
            disp(['Not creating path structure since there are files or directories in ' fetpth]);
        else
            eval(['!(cd ' fetpth '; (cd ' par.speechpth '; find * -type d) |xargs mkdir -p)']);
        end;
        eval(['!(cd ' par.speechpth '; find * -type f -name "*.' frm '") > tmp1.scp;']);    % find all speech files
        if ~isempty(grepstr)
            eval(['!egrep ' grepstr ' tmp1.scp > tmp2.scp']);                               % only keep those that match grepstr
            eval(['!rm -f tmp1.scp']);
        else
            eval(['!mv tmp1.scp tmp2.scp']);                                                % keep all files if grepstr is empty
        end;
        eval(['!sed "s:^:' fetpth '/:" tmp2.scp > tmp3.scp']);
        eval(['!rm -f tmp2.scp']);
        eval(['!sed "s:' frm ':' fetext ':" tmp3.scp > ' fetlistname]);
        eval(['!rm -f tmp3.scp']);
    end
end

