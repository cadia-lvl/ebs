function mt=fftresm(m,n,p)
% create matrix to resample FFT in order to change the frequency resolution
%
%  Inputs: m        size of input FFT
%          n        size of output FFT
%          p        odd order of the filter generating the intermedate signal [5]
%  Output: mt(m,n)  sparse transformation matrix whose columns sum to 1
%
% Conceptually, the routine first converts the input spectrum to an
% intermediate piecewise linear spectrum and then averages this over
% the output frequency bin widths (0:n-1)/n+-0.5/n to get the output
% spectrum. The intermediate spectrum is chosen such that if n=m, the
% output matrix will equal the identity.
%
% Note: it is possible that we can make a better FIR pre-filter of a given order using, e.g. Remezexchange algorithm. However, this would depend on m and would therefore be much slower.
%
persistent m0 n0 mt0 p0 a0
if nargin<3
    p=5;
else
    p=p+1-mod(p,2);                             % round p up to the next higher odd number
end
if isempty(m0) || m~=m0 || n~=n0 || p~=p0
    p2=0.5*(p-1);                                   % number of elements in first half of pre-filter excluding central element
    if m>=n                                                         % we are decreasing the spectral resolution. Input segment i goes between input i to i+1 (where DC=1)
        klo=repmat((0:m-1)/m,1,1).*repmat(n,1,m)+0.5;               % input node positions in shifted fractional output units (where DC output is 0.5)
        kkmid=ceil(klo);                                            % input segment between i=1:m contributes to outputs kkmid and kkmid+1 (mod n)
        fmid=(kkmid-0.5)*m/n-(0:m-1);                               % distance from each input node to next output integration boundary
        fhi1=min(fmid,1);                                           % fraction of input segment i=1:m that contribiutes to output kkmid
        wlo2=0.5*fhi1.^2;                                           % contribution of input i+1 (mod m) in segement i to the output kkmid
        wlo1= fhi1-wlo2;                                            % contribution of input i to the output kkmid
        tt=sparse([1:m 1:m 2:m 1 2:m 1],mod([kkmid-1 kkmid kkmid-1 kkmid],n)+1,(n/m)*[wlo1 0.5-wlo1 wlo2 0.5-wlo2],m,n); % <=3 non-zero entries per row
    else                                                                        % we are increasing the spectral resolution
        % each output bin is formed from at most two input bins
        jlo=repmat(((0:n-1)-0.5)/n,1,1).*repmat(m,1,n);                         % low edge of each output's integration range in fractional input units (DC = 0)
        jjmid=ceil(jlo);                                                        % round up to next higher input mid-bin frequency
        kwid=m/n;                                                               % output bin width in fractional input units
        khi=jlo+kwid;                                                           % high edge of integration range in fractional input units
        flo1=jlo-jjmid;                                                         % lower integration limit (always -ve)
        fhi=khi-jjmid;                                                          % upper integration limit; maybe -ve (only one input bin contributes) or +ve (two input bins contribute)
        fhi1=min(fhi,0);                                                        % upper limit of first integration (never +ve)
        fhi2=max(fhi,0);                                                        % upper limit of second integration (never -ve). Note that lower limit is always 0
        % note fhi2+fhi1-flo1 = repmat(kwid,1,n)
        wfact=0.5/kwid;                                                         % weight normalization factor
        wlo=wfact*(flo1.^2-fhi1.^2);                                            % weight for jjlo input
        whi=wfact*fhi2.^2;                                                      % weight for jjhi input
        tt=sparse(mod([jjmid-1 jjmid jjmid+1],m)+1,repmat(1:n,1,3),[wlo 1-wlo-whi whi],m,n); % <=3 entries per column
    end
    if isempty(p0) || p~=p0
        % If all frames are the same length, and flout=flin, then the output signal is the intermediate signal filtered by [1 6 1]/8.
        % To calculate the intermediate signal, we therefore would like to apply a filter of order p  that, when convolved with  [1 6 1]/8
        % equals an impulse. This is not possible, so we instead seek a least-squares solution and apply a constraint to force unity DC gain.
        % Since the filter is assumed symmetric, we only calculate its first half.
        if p<3
            a=1;                                        % if p=1, the pre-filter is just a delta function
        else
            fm=0.125*sqrt(2)*toeplitz([1; 6; 1; zeros(p2-1,1)],[1 zeros(1,p2)]); % upper half of forward filter matrix of size p2+2,p2+1) (falsely assuming constant-length frames)
            fm(end,end-1:end)=[0.25 0.75];
            b=[zeros(p2+1,1);1];
            k=2-b(2:end);                               % constraint vector = [2; ... 2; 1] (sum of coefficients = 1)
            ahl=[fm'*fm k; k' 0]\[fm'*b; 1];            % calculate first half of filter and the Lagrange multiplier
            a=[ahl(1:end-1);ahl(end-2:-1:1)];           % symmetric FIR pre-filter of order p
        end
        p0=p;
        a0=a;
    else
        a=a0;
    end
    tpr=repmat((1:m)',1,p);                         % pre-filter matrix row index
    mt=sparse(tpr,mod(tpr+repmat(-p2:p2,m,1)-1,m)+1,repmat(a',m,1),m,m)*tt; % transformation matrix
    m0=m;
    n0=n;
    mt0=mt;
else
    mt=mt0;
end
if ~nargout
    imagesc([0 m-1]/m, [0 n-1]/n, mt');
    axis 'xy';
    colorbar;
    hold on;
    plot([0 1],[0 1],'-k+'); % plot the diagonal
    hold off
    xlabel(sprintf('Input Freq (%d)',m));
    ylabel(sprintf('Output Freq (%d)',n));
    cblabel('Weight');
end
