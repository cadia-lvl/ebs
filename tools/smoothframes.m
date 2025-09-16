function feout=smoothframes(s,fein,p)
% Smooth framelengths across unvoiced intervals
%
%  Inputs:      s(n,1)  Speech signal
%            fein(1,m1) last sample number in each frame
%               p       optional structure giving procesing options:
%                           p.smoothalg = {'none','lin','quadlin','quadlog'} ['none']
%                           p.voithresh = cross-correlation threshold for voiced frames [0.6]
%                           p.voilenrat = length ratio threshold for voiced frames [1.1]
%
% Outputs:  feout(1,m2) revised last sample number in each frame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   algorithm Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent q0 tqlin
if isempty(q0)
    q0.smoothalg='none';        % smoothing algorithm {'none','lin','quadlin','quadlog'}
    q0.voithresh=0.6;           % cross-correlation threshold for voiced frames
    q0.voilenrat=1.1;           % maximum length ratio of successive voiced frames
    tqlin=inv([0 0 1 0; 1 1 1 0; 0 1 0 -1; 2 1 0 -1]); % save transformation matfsinrix for 'quadlin' smoothing
end
if nargin<3
    q=q0;
else
    q=v_paramsetch(q0,p);             % update parameters from p input
end
feout=fein; % by default just copy input framing to output
if ~strcmp(q.smoothalg,'none')   % if smoothing required
    nfin=length(fein); % number of frames
    flin=filter([1 -1],1,fein); % frame lengths
    fsin=fein-flin+1; % first sample of each frame
    minflin=min(flin(2:end),flin(1:end-1)); % minimum of consecutive frame lengths
    xcrin=zeros(1,nfin-1); % space for inter-frame correlation coefficients; xcrin(i) is correlation between frames i and i+1
    for i=1:nfin-1 % could possibly do this without loops with some fancy indexing
        s1=s(fsin(i):fsin(i)+minflin(i)-1); % extract minflin samples from first frame
        s2=s(fsin(i+1):fsin(i+1)+minflin(i)-1);  % extract  minflin samples from second frame
        s1=s1-mean(s1);
        s2=s2-mean(s1);
        xcrin(i)=(s1(:)'*s2(:))/sqrt(s1(:)'*s1(:)*s2(:)'*s2(:)); % calculate correlation coefficient
    end
    %%%%%%%%%%% Select valid frames %%%%%%%%%%%%%
    % for now we just use the adjacent-frame cross-correlation but we could also use frame length ratio or require multiple consecutive valid frames
    lrat=flin(2:end)./flin(1:end-1); % length ratio of adjacent frames
    xcrv=(xcrin>=q.voithresh) & (lrat<=q.voilenrat) & (lrat>=1/q.voilenrat); % find frames with good x-correlation an a good length ratio
    frameok=[xcrv true] & [true xcrv]; % mask for valid frames
    segstart=find(~frameok & [true frameok(1:end-1)]); % first of string of bad frames
    segend=find(~frameok & [frameok(2:end) true]); % last of string of bad frames
    if ~isempty(segstart) && (segstart(1)>1 || segend(1)<nfin) % check if anything to do
        feout=[];                           % initialize feout
        jb=1;                               % dummy last frame of previous bad string
        for i=1:length(segstart)            % loop for each string of consecutive bad frames
            ja=segstart(i);                 % first frame in string
            feout=[feout fein(jb:ja-1)];    % copy previous string of good frame boundaries to the output
            jb=segend(i);                   % last frame in string
            ta=fsin(ja)-1;                  % end of previous good frame
            tgap=fein(jb)-ta;               % length of gap (samples)
            if ja==1                        % if string starts with first frame
                fa=tgap/flin(jb+1);         % assume pitch of previous good frame equals that of next good frame
            else
                fa=tgap/flin(ja-1);         % scaled pitch of previous good frame (in cycles/gap)
            end
            if jb==nfin                     % if string ends with final frame
                fb=fa;                      % assume pitch of next good frame equals that of previous good frame
            else
                fb=tgap/flin(jb+1);         % scaled pitch of next good frame (in cycles/gap)
            end
            k=round(0.5*(fa+fb));           % number of periods in gap is k ( possibly <1 )
            if k>1                          % we need to insert some extra frame boundaries
                switch q.smoothalg
                    case 'lin'                              %%%%% 'lin': insert equal length frames %%%%%%
                        tadd=ta+round((1:k-1)*tgap/k);      % equally spaced pseudo-epochs
                    case 'quadlin'                          %%%%% 'quadlin': pitch increases linearly with equal pitch error, p(4), at each end %%%%
                        p=tqlin*[0 k fa fb]';               % frame-phase, q, = p(1)*t^2+p(2)*t+p(3) for 0<t<1, so instantaneous pitch = dq/dt=2*p(1)*t+p(2)
                        kk=1:k-1;
                        tadd=ta+round(tgap*2*(p(3)-kk)./(-sqrt(p(2)^2-4*p(1)*(p(3)-kk))-p(2)));        % solve for frame-phase = integers kk
                    case 'quadlog'                          %%%%% 'quadlin': pitch increases linearly with equal log-pitch error at each end %%%%
                        p=[0 0 1 0; 1 1 1 0; 0 1 0 -fa; 2 1 0 -fb]\[0 k 0 0]';                  % frame-phase = p(1)*t^2+p(2)*t+p(3)
                        kk=1:k-1;
                        tadd=ta+round(tgap*2*(p(3)-kk)./(-sqrt(p(2)^2-4*p(1)*(p(3)-kk))-p(2)));        % solve for frame-phase = integers kk
                end
            else                                            % no additional epochs to insert
                tadd=[];
            end
            feout=[feout tadd]; % insert extra frames
        end
    end
end
if ~nargout
    % if strcmp(q.smoothalg,'none')                   % if no smoothing, we must calculate xcrin for plotting
    %     nfin=length(fein);                          % number of input frames
    %     flin=filter([1 -1],1,fein);                 % frame lengths (samples)
    %     fsin=fein-flin+1;                           % first sample of each frame
    %     minflin=min(flin(2:end),flin(1:end-1));     % minimum of consecutive frame lengths
    %     xcrin=zeros(1,nfin-1);                      % space for inter-frame correlation coefficients; xcrin(i) is correlation between frames i and i+1
    %     for i=1:nfin-1                              % could possibly do this without loops with some fancy indexing
    %         s1=s(fsin(i):fsin(i)+minflin(i)-1);     % extract minflin samples from first frame
    %         s2=s(fsin(i+1):fsin(i+1)+minflin(i)-1); % extract  minflin samples from second frame
    %         s1=s1-mean(s1);                         % subtract means
    %         s2=s2-mean(s1);
    %         xcrin(i)=(s1(:)'*s2(:))/sqrt(s1(:)'*s1(:)*s2(:)'*s2(:)); % calculate correlation coefficient for adjacent input frames
    %     end
    % end
    nfout=length(feout);                            % number of output frames
    flout=filter([1 -1],1,feout);                   % frame lengths
    fsout=feout-flout+1;                            % first sample of each frame
    minflout=min(flout(2:end),flout(1:end-1));      % minimum of consecutive frame lengths
    xcrout=zeros(1,nfout-1);                        % space for inter-frame correlation coefficients; xcrin(i) is correlation between frames i and i+1
    for i=1:nfout-1                                 % could possibly do this without loops with some fancy indexing
        s1=s(fsout(i):fsout(i)+minflout(i)-1);      % extract minflin samples from first frame
        s2=s(fsout(i+1):fsout(i+1)+minflout(i)-1);  % extract  minflin samples from second frame
        s1=s1-mean(s1);
        s2=s2-mean(s1);
        xcrout(i)=(s1(:)'*s2(:))/sqrt(s1(:)'*s1(:)*s2(:)'*s2(:)); % calculate correlation coefficient for adjacent output frames
    end
    axlinkt=[];
    %%%%%%%%
    if strcmp(q.smoothalg,'none')                   % if no smoothing, the input and output plots will be the same
        subplot(3,1,3);                                 %%%%% Plot Adjacent frame length ratio
        semilogy(feout(1:end-1),flout(2:end)./flout(1:end-1),'-r');
        v_axisenlarge([-1 -1.05]);
        xlim=get(gca,'xlim');
        hold on;
        plot(xlim,[1 1]*q.voilenrat,'--k',xlim,[1 1]/q.voilenrat,'--k');
        hold off
        ssi=yticksi;
        ylabel('Adjacent Frame Ratio');
        axlinkt=[axlinkt gca]; % link time axis
        %%%%%%%%
        subplot(3,1,2);                                 %%%%% Plot Adjacent frame correlations
        plot(feout(1:end-1),xcrout,'-r');
        v_axisenlarge([-1 -1.05]);
        xlim=get(gca,'xlim');
                hold on;
        plot(xlim,q.voithresh*[1 1],'--k');
        hold off
        ylabel('Adjacent Frame Corr');
        axlinkt=[axlinkt gca]; % link time axis
        %%%%%%%%
        subplot(3,1,1);                                 %%%%% Plot Frame lengths
        ttout=[fsout;feout];                            % start and end of output frames (samples)
        llout=repmat(flout,2,1);
        plot(ttout(:),llout(:),'-r');
        v_axisenlarge([-1 -1.05]);
        ylabel('Frame Len (samp)');
        axlinkt=[axlinkt gca]; % link time axis
    else
        subplot(3,1,3);                                 %%%%% Plot Adjacent frame length ratio
        semilogy(fein(1:end-1),flin(2:end)./flin(1:end-1),'-b',feout(1:end-1),flout(2:end)./flout(1:end-1),'-r');
        v_axisenlarge([-1 -1.05]);
        xlim=get(gca,'xlim');
        hold on;
        plot(xlim,[1 1]*q.voilenrat,'--k',xlim,[1 1]/q.voilenrat,'--k');
        hold off
        legend('Original','Smoothed','Location','Best');
        ssi=yticksi;
        ylabel('Adjacent Frame Ratio');
        axlinkt=[axlinkt gca]; % link time axis
        %%%%%%%%
        subplot(3,1,2);                                 %%%%% Plot Adjacent frame correlations
        plot(fein(1:end-1),xcrin,'-b',feout(1:end-1),xcrout,'-r');
        v_axisenlarge([-1 -1.05]);
        xlim=get(gca,'xlim');
        hold on;
        plot(xlim,q.voithresh*[1 1],'--k');
        hold off
        legend('Original','Smoothed','Location','Best');
        ylabel('Adjacent Frame Corr');
        axlinkt=[axlinkt gca]; % link time axis
        %%%%%%%%
        subplot(3,1,1);                                 %%%%% Plot Frame lengths
        ttin=[fsin;fein];                               % start and end of input frames (samples)
        llin=repmat(flin,2,1);                          % input frame lengths (samples)
        ttout=[fsout;feout];                            % start and end of output frames (samples)
        llout=repmat(flout,2,1);
        plot(ttin(:),llin(:),'-b',ttout(:),llout(:),'-r');
        v_axisenlarge([-1 -1.05]);
        legend('Original','Smoothed','Location','Best');
        ylabel('Frame Len (samp)');
        axlinkt=[axlinkt gca]; % link time axis
    end
    title(['Epoch Smoothing (' q.smoothalg ')']);
    linkaxes(axlinkt,'x');
end