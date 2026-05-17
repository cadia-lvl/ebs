function plotstftemeta(mode,meta,fs,scale,offset)
%  Inputs: mode         char string with options 
%          meta(:,10)   metadata
%          fs           sample frequency in Hz [1]
%          scale        multiply plot by scale factor and then add offset (-x = scale [0,x])
%          offset       multiply plot by scale factor and then add offset
%
% mode options: 'l'     frame length
%               'f'     Inverse frame length (approx pitch) [NOT IMPLEMENTED]
%               'n'     DFT length
%               'o'     offset [NOT IMPLEMENTED]
%               's'     scale factor [NOT IMPLEMENTED]
%               'g'     group delay
%               'e'     EWGD
%               'p'     EWPD
%               'P'     DC Phase [NOT IMPLEMENTED]
%               'G'     DC group delay
%               'w'     plot frames as waveform ('W' to invert or 'wW' to centre)
%
%               't'     plot values in seconds (or 1/Hz) rather than samples
%               'T'     plot values in milliseconds (or 1/kHz) rather than samples
%               'v'     make frame boundaries midway between samples
%               'z'     first sample is at t=0 rather than t=1/fs
%             
%               ';'     chars follwing ';' are given in the plot command [';-k']
%
% meta(*,:)=[first-sample, frame-length, dft-length, offset, scale-factor, group-delay (samples), EWGD (samples), EWPD (samples), DC phase (rad), DC group delay (samples)]
%
if nargin<3 || isempty(fs)
    fs=1;
end
if nargin<4 || isempty(scale)
    scale=1;
end
if nargin<5 || isempty(offset)
    offset=0;
end
mode=char(mode);                % force mode to be char array
co=strfind(mode,';');
if isempty(co)
    modeg=mode;
    modep='';
else
    modeg=mode(1:co(1)-1);      % mode string for graph selection
    modep=mode(co(1)+1:end);    % mode string for plotting
end
if isempty(modep)
    modep='-';
end
nframe=size(meta,1); % number of frames
if ~isempty(modeg) || nframe==0 % nothing to plot if modeg is empty or no frames
    vopt=0.5*any(modeg=='v');
    fstart=meta(:,1)-vopt;
    fend=fstart+meta(:,2)-(1-2*vopt);
    tplotopt='lngepG'; % time plot options
    ntplotopt=length(tplotopt);
    tplotleg=["Frame Len" "DFT Len" "Group Del" "EWGD" "EWPD" "DC Grp Del"];
    tplotcol=[2 3 6 7 8 10];
    ntplot=sum(reshape(repmat(modeg',1,ntplotopt)==repmat(tplotopt,length(modeg),1),[],1)); % number of time-value quantities to plot
    if ntplot>0 % y axis is samples or seconds
        leg=repmat(" ",1,ntplot); % legend array
        itplot=0;
        taxis=(reshape([fstart'; fend'],[],1)-any(modeg=='z'))/fs;
        v=zeros(2*nframe,ntplot); % space for plot values
        for i=1:ntplotopt
            if any(modeg==tplotopt(i))
                itplot=itplot+1;
                leg(itplot)=tplotleg(i);
                v(:,itplot)=reshape(repmat(meta(:,tplotcol(i))',2,1),[],1);
            end
        end
        if scale>0
            v=v*(scale*(1+999*any(modeg=='T'))/(1+(fs-1)*any(lower(modeg)=='t')))+offset;
        else % normalize scaling
            minv=min(v(:));
            maxv=max(v(:));
            v=(v-minv)*(scale/(minv-maxv-(maxv==minv)))+offset;
        end
        plot(taxis,v,modep)
        v_axisenlarge([-1 -1.05]);
        if fs==1
            xlabel('Samples');
        else
            xlabel('Time (s)');
        end
        if scale==1 && offset==0
            if any(modeg=='T')
                tunit='ms';
            elseif any(modeg=='t')
                tunit='s';
            else
                tunit='samples';
            end
        else
            tunit='';
        end
        if ntplot>1
            legend(leg,'location','best');
            if isempty(tunit)
                ylabel('Length');
            else
                ylabel(['Length (' tunit ')']);
            end
        else
            if isempty(tunit)
                ylabel(leg(1));
            else
                ylabel([char(leg(1)) ' (' tunit ')']);
            end
        end
    elseif any(lower(modeg)=='w')
        % send=sum(meta(:,1:2))-1; % last sample in each frame
        % levs=zeros(nframe,1); % assigned level of each frame
        % levend=zeros(nframe,1); % end of last frame on each level
        % nxtlev=1; % next available level
        % for ifr=1:nframe % assign each frame to a level
        %     ilev=find(levend(1:nxtlev)<meta(ifr,1),1)
        %     nxtlev=nxtlev+(levend(ilev)==0); % update the numbewr of levels
        %     levs(ifr)=ilev; % assign this frame to ilev
        %     levend(ilev)=send(ifr);
        % end
    end
end

