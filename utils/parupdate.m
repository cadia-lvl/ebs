function [p,descrip]=parupdate(pin,config,posfields,format)
%PARUPDATE update parameter structure to try different configurations
%
% Usage: par0=struct('a',1,'b',2,'c','three');                      % default parameter structure
%        ff={'c' 'a'};                                              % cell array listing the position-dependent parameters within configpars structure
%        config={                                                   % one row per trial giving the parameter values
%                 {'1st' 11};                                       % ... parameters for trial 1
%                 {'2nd' 1 'd' 7};                                  % ... parameters for trial 2 includes new parameter: 'd'
%               }; 
%        [par,desc]=parupdate(par0,config{1},ff);                   % desc='a=11, c=1st'
%        [par,desc]=parupdate(par0,config{1},ff,'%c%: %a%,%b%');    % desc='1st: 11,2'
%        [par,desc]=parupdate(par0,config{2},ff);                   % desc='c=2nd, d=7' (only includes fields changed from default values)
%
%  Inputs:  pin         inpout parameter structure with default values
%           config      value of fixed fields followed by name-value pairs
%           posfields   list of position-dependent fields in config input
%           format      format of descrip output [default: '$c&=%&%$, $']
%                         (a) %field% interpolates field value
%                         (b) $x...$...$ loops over fields selected by hex digit x which
%                             is the sum of: 1=default, 2=unchanged, 4=changed, 8=added. The first segment
%                             will be included for each field, whereas the second segment (if present) for
%                             each field *except* the last.
%                         (c) & or %&% interpolates the name or value of the current loop's field.
%                         (d) \%, \&, \$, \\ escape special chars.
%
% Outputs:  p           output structure with updated field values
%           descrip     text describing non-defaut paramerter values
%
p=pin;
nfixed=0;
if nargin>2
    nfixed=length(posfields);
    for i=1:nfixed
        p.(posfields{i})=config{i};
    end
end
if nargin>1 % some fields need updating
    for i=nfixed+1:2:length(config)                             % loop through list of updatable parameters
        p.(config{i})=config{i+1};                            % update name-value pairs
    end
end
if nargout~=1                                                  % descriptive text required
    descrip='';
    fnames=fieldnames(p);                                     % list the field names in the updated structure
    nfnames=length(fnames);                                     % number of field names
    fnamstat=zeros(nfnames,1);                                  % status of fields
    fnamstat(:)=1;
    for i=1:nfixed
        fnm=posfields{i};                                       % field name
        msk=strcmp(fnm,fnames);                                 % mask identifiying this field
        if isfield(pin,fnm)                                   % updated field
            fnamstat(msk)=2*(2-isequal(pin.(fnm),p.(fnm)));   % 2=same value, 4=new value
        else                                                    % newly created field
            fnamstat(msk)=8;                                  % 8=new field
        end
    end
    if nargin>1 && length(config)>nfixed+1                     % some fields were updated
        for i=nfixed+1:2:length(config)                         % loop through list of updatable parameters
            fnm=config{i}; % field name
            msk=strcmp(fnm,fnames);                                 % mask identifiying this field
            if isfield(pin,fnm)                                   % updated field
                fnamstat(msk)=2*(2-isequal(pin.(fnm),p.(fnm)));   % 2=same value, 4=new value
            else                                                    % newly created field
                fnamstat(msk)=8;                                  % 8=new field
            end
        end
    end
    if nargin<4
        format={'$c&=%&%$, $'}; % default to listing updated fields
    end
    descrip=repmat(' ',1,1000); % initialize description to all spaces
    j=0; % initialize index into descrip
    i=0; % index into format
    state=0; % initialize state machine state
    nformat=length(format);
    while i<nformat
        i=i+1;
        fmi=format(i);
        if fmi=='\' && i<nformat % escape next character
            j=j+1;
            i=i+1;
            descrip(j)=format(i); % copy character to output
        else
            switch state
                case 0 % normal text
                    switch fmi
                        case '%' % interpolate field value
                            state=1;
                            fieldname=repmat(' ',1,100);            % initialize field name
                            k=0;                                    % and index
                        case '$'                                    % loop of fields
                            idoll=find(format(i:end)=='$')+i-1;     % find loop boundaries
                            if length(idoll)>1
                                idoll=idoll([1 2 min(3,length(idoll))]); % positions of three $ that define loop (last two might be coincident)
                                state=2;
                            end
                        otherwise
                            j=j+1;
                            descrip(j)=fmi; % copy character to output
                    end

                case 1 % read expression to interpolate
                    switch fmi
                        case '%' % end of field name
                            fieldid=find(strcmp(fieldname(1:k),fnames));
                            if k>0 && ~isempty(fieldid)
                                val=char(strtrim(formattedDisplayText(p.(fieldname(1:k))))); % value of field
                                descrip(j+1:j+length(val))=val;
                                j=j+length(val);
                                fnamstat(fieldid)=0; % remove field from possible future loop
                                state=0;
                            end
                        otherwise % partway through field name
                            k=k+1;
                            fieldname(k)=fmi;
                    end

                case 2          % read new loop specification
                    fmi=lower(fmi); % convert hex digit to lower case
                    loopspec=fmi-48-(fmi>96)*39; % bit mask: 1=default, 2=unchanged, 4=changed, 8=added
                    fieldlist=find(bitand(loopspec,fnamstat)>0); % list of fields to display
                    nfield=length(fieldlist);
                    if nfield>0
                        ifield=1; % initialize field index
                        state=3; % within a loop
                    else % no loop entries to process
                        i=idoll(3); % skip to end of loop
                        state=0;
                    end
                case 3          % normal text within a loop
                    switch fmi
                        case '%' % interpolate field value
                            state=4;
                            fieldname=repmat(' ',1,100); % initialize field name
                            k=0; % and index
                        case '$' % loop of fields
                            if i==idoll(3) % end of loop
                                if ifield<nfield
                                    ifield=ifield+1;
                                    i=idoll(1)+1; % reset loop pointer skipping field specifier
                                else
                                    state=0; % end of loop. Continue.
                                end
                            else
                                if ifield==nfield % if this is the final time through the loop
                                    i=idoll(3); % skip to end of loop
                                    state=0;
                                end
                            end
                        case '&' % interpolate current field name
                            val=fnames{fieldlist(ifield)}; % value of field
                            descrip(j+1:j+length(val))=val;
                            j=j+length(val);
                        otherwise
                            j=j+1;
                            descrip(j)=fmi; % copy character to output
                    end

                case 4 % read field name to interpolate within a loop
                    switch fmi
                        case '%' % end of field name
                            if strcmp(fieldname(1:k), '&')
                                fieldid=fieldlist(ifield); % get current field ID
                            else
                                fieldid=find(strcmp(fieldname(1:k),fnames));
                            end
                            if k>0 && ~isempty(fieldid)
                                val=char(strtrim(formattedDisplayText(p.(fnames{fieldid}),'NumericFormat','shortG'))); % value of field
                                descrip(j+1:j+length(val))=val;
                                j=j+length(val);
                                fnamstat(fieldid)=0; % remove field from possible future loop
                                state=3; % return to loop processing
                            end
                        otherwise % partway through field name
                            k=k+1;
                            fieldname(k)=fmi;
                    end
            end
        end
    end
    descrip=descrip(1:j);
end
if ~nargout
    fprintf(1,'%s\n',descrip);
end