function stfteplot(stfte,meta,fs,dbr)
% plot stfte as a spectrogram
%  Inputs:  stfte(nframe,maxbin)    Complex stfte data
%           meta(:,6)               meta(i,1:3)=[first-sample frame-length DFT-length] in samples for frame i
%           dbr                     range in dB either [min max] or [range] or a linear scale: 0=abs value, -1=real part, -2=imag part [default=40]
if nargin<4 || isempty(dbr)
    dbr=40;
end
if nargin<3 || isempty(fs)
    fs=1;
end
dblin=numel(dbr)==1 && dbr<=0; % for a linear scale:  0=abs value, -1=real part, -2=imag part
[nframe,maxbin]=size(stfte);
frc=meta(:,1:2)*[1;.5]-0.5;                                                     % frame centres
fre=0.5*[[3 -1]*frc(1:2); frc(1:end-1)+frc(2:end); [-1 3]*frc(end-1:end)];      % frame edges (needs >=2 frames)
nmap=size(colormap,1);                                                          % number of colormap entries
maxbinp=1+floor(maxbin/2);                                                      % max number of positive frequencies
patylo=max(repmat(-1:2:maxbin-1,nframe,1)./(2*repmat(meta(:,3),1,maxbinp)),0);  % lower edge of each pach
patyhi=min(repmat(1:2:maxbin+1,nframe,1)./(2*repmat(meta(:,3),1,maxbinp)),0.5); % upper edge of each pach
patxlo=repmat(fre(1:nframe),1,maxbinp);                                         % left edge of each pach
patxhi=repmat(fre(2:nframe+1),1,maxbinp);                                       % right edge of each pach
patxv=repmat(1:maxbinp,nframe,1)<=1+floor(meta(:,3)/2);                         % Valid datapoints corresponding to +ve frequencies
patx=[patxlo(patxv)'; patxhi(patxv)'; patxhi(patxv)'; patxlo(patxv)'];          % x coordinate of patch corners
paty=[patylo(patxv)'; patylo(patxv)'; patyhi(patxv)'; patyhi(patxv)'];          % y coordinate of patch corners
if dblin                                                                        % if linear value scale
    switch dbr
        case -1
            patcm=real(stfte(:,1:maxbinp));                                              % patch value as real paRT
                      clim=[min(patcm(:)) max(patcm(:))];                                                     % range is 0 to max
            if isreal(stfte)
                clab="Value";
            else
                clab="Real Part";
            end
        case -2
            patcm=imag(stfte(:,1:maxbinp));                                              % patch value as imag part
                   clim=[min(patcm(:)) max(patcm(:))];                                                     % range is 0 to max
            clab="Imag Part";
        otherwise
            patcm=abs(stfte(:,1:maxbinp));                                              % patch value as absolute value
            clim=[0 max(patcm(:))];                                                     % range is 0 to max
            clab="Value";                                                               % colorbar label
    end
else
    patcm=db(stfte(:,1:maxbinp));                                               % patch value in dB
    if numel(dbr)>1                                                             % lower and upper dB values are specified
        clim=[dbr(1) dbr(2)];
    else                                                                        % dB range is specified
        clim=[-dbr(1) 0]+max(patcm(:));
    end
    clab="Value (dB)";
end
cla;                                                                            % clear figure: should maybe check first if "hold" is on
patch(patx/fs,paty*fs,patcm(patxv)','EdgeColor','none');
cbh=colorbar;
set(gca,'CLim',clim,'xlim',[fre(1) fre(end)]/fs,'ylim',[0 0.5*fs]);
colormap("parula");
cblabel(clab);
if fs==1
    xlabel('Samples');
    ylabel('Frequency\div{}f_s');
else
    xlabel(['Time (' v_xticksi 's)']);
    ylabel(['Frequency (' v_yticksi 'Hz)']);
end