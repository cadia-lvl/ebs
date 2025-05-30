function p=gettimitpath
%GETTIMITPATH get path of timit subfolder of timit database on the current computer
persistent q
if ~isempty(q)
    p=q;                        % use previously cached path if possible
else
    switch computername
        case 'mac.home'
            p='/Users/jg/Data/TIMIT/timit/';
        case 'ee-dmb4'
            p='D:/OneDrive - Imperial College London/work/data/speech/timit/timit/';    % path to timit sub-folder of timit CD
        case 'elitedesk'
            p='D:/OneDrive - Imperial College London/work/data/speech/timit/timit/';    % path to timit sub-folder of timit CD
        otherwise
            error(sprintf('Add your computername (''%s'') into the switch statement in %s.m',computername,mfilename));
    end
    q=p;                        % save path in cache for future calls
end