function [framex,costgaindb,finaloff]=gs_frameadj(framek,s,fs,maxoff,maxlen)
%GS_FRAMEADJ adjust frame boundaries to minimize squared boundary discontinuity
%
%  Inputs: framek(1,n)   Last sample numbers in each of n frames. Hence framek(end)=length(s).
%               s(:,1)   Speech signal
%                   fs   Sample frequency
%               maxoff   Maximum permitted frame-boundary adjustment in samples
%               maxlen   Maximum frame length in samples [Inf]
%
% Outputs: framex(1,n)  Last sample numbers in each of n frames. Hence framex(end)=length(s).
%           costgaindb  Improvement in mean squared boundary discontinuity in dB
%        finaloff(1,n)  Number of samples by which each boundary has been changed. Hence finaloff(end)=0.
%
nframe=length(framek);
s=s(:);                         % force s to be a column vector
if nargin<5
    maxlen=Inf;                 % no maximum length
end
%
% adjust frames slightly to minimize transients
%
finaloff=zeros(1,nframe);
costinit=sum((s(framek(2:end))-s(framek(1:end-1))).^2); % sum of squared boundary discontinuity
framex=framek; % copy original frame end samples
if maxoff>0
    offs=-maxoff:maxoff;                                        % list of possible boundary offsets
    noffs=length(offs);                                         % number of possible offsets
    prev=zeros(nframe,noffs);                                   % best previous offset | current offset
    cost=[0 ones(1,noffs-1)];                                   % cost of best path given current boundary offset
    soff=repmat(s(framek(1)),1,noffs);                          % initialise all paths to final sample of frame 1 (as surrogate for frame 0)
    woff=ones(1,noffs);                                         % repeated 1's used for indexing
    lenmat=offs(woff,:)-offs(woff,:)';
    for i=2:nframe-1                                            % Dynamic programming loop: current boundary is last sample of frame i
        oldsoff=soff';                                          % final sample of previous frame | previous offset
        soff=s(framek(i)+offs)';                                 % final sample of current frame | current offset
        oldcost=cost';                                          % cost of best path | previous offset
        % cost terms: previous cost + squared difference between final samples + Inf penalty for exceeding maxlen
        [cost,prev(i,:)]=min(oldcost(:,woff)+(soff(woff,:)-oldsoff(:,woff)).^2+(1./(lenmat<=maxlen+framek(i-1)-framek(i))-1),[],1);
    end
    [dum,pr]=min(cost);
    for i=nframe-1:-1:2                                         % do dynamic programming traceback
        finaloff(i)=offs(pr);
        framex(i)=framek(i)+offs(pr);                           % adjust frame boundary
        pr=prev(i,pr);
    end
end
costfinal=sum((s(framex(2:end))-s(framex(1:end-1))).^2);        % sum of squared boundary discontinuity
costgaindb=db(costinit/costfinal)/2;                            % cost improvement in dB
