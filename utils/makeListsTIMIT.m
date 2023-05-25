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

if ispc
    error('Not implemented yet for Windows');
end;

if nargin<2
    grepstr='';
end;

if nargin<2
    grepstr='';
end;

%Short-hands
opt=par.listopt;
frm=par.frm;
splistname=par.splistname;
fetlistname=par.fetlistname;
phlistname=par.phlistname;
fetext=par.fetext;


if any(opt=='s')  %Make splist
    eval(['!(cd ' par.speechpth '; find * -type f -name "*.' frm '") > tmp1.scp;']);
    if ~isempty(grepstr)
        eval(['!egrep "' grepstr '" tmp1.scp > tmp2.scp']);
        eval(['!rm -f tmp1.scp']);
    else
        eval(['!mv tmp1.scp tmp2.scp']);
    end;
    eval(['!sed "s:^:' par.speechpth '/:" tmp2.scp > ' splistname]);
    eval(['!rm -f tmp2.scp']);   
end;

if any(opt=='p')
    eval(['!sed "s:' frm ':PHN:" ' splistname ' > ' phlistname]);
end;
    
if any(opt=='f')
    dd=dir(par.fetpth);
    if length(dd)>2
        disp(['Not creating path structure since there are files or directories in ' par.fetpth]);
    else
        eval(['!(cd ' par.fetpth '; (cd ' par.speechpth '; find * -type d) |xargs mkdir -p)']);
    end;
    eval(['!(cd ' par.speechpth '; find * -type f -name "*.' frm '") > tmp1.scp;']);
    if ~isempty(grepstr)
        eval(['!egrep ' grepstr ' tmp1.scp > tmp2.scp']);
        eval(['!rm -f tmp1.scp']);
    else
        eval(['!mv tmp1.scp tmp2.scp']);
    end;
    eval(['!sed "s:^:' par.fetpth '/:" tmp2.scp > tmp3.scp']);
    eval(['!rm -f tmp2.scp']);
    eval(['!sed "s:' frm ':' fetext ':" tmp3.scp > ' fetlistname]);
    eval(['!rm -f tmp3.scp']);
end;

