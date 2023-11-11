function inm=lininterps(x,y,md)
%  Inputs: x(n)   sample points of input function
%          y(m)  sample points of output function
%          md       character string containing a combination of
%                     'z','Z'   Set to zero at low/high ends [default]
%                     'r','R'   replicate first/last x value at low/high ends
%                     'e','E'   linearly extrapolate at low/high ends
%
% Outputs: inm(m,n)   sparse interpolation matrix: vy=inm*vx
%
n=length(x); % number of input points
m=length(y); % number of output points
mn=m+n;
if nargin<3 || isempty(m)
    ws=0; % default end processing
    we=0;
else
    ws=min(any(repmat('re',length(md),1)==repmat(md',1,2),1)*(1:2)',2); % start processing: 0=zero, 1=replicate, 2=extrapolate
    we=min(any(repmat('RE',length(md),1)==repmat(md',1,2),1)*(1:2)',2); % end processing: 0=zero, 1=replicate, 2=extrapolate
end
[fint,ifint]=sort([y(:)' x(:)']);           % intertwine input and output knots
jfint=zeros(1,mn);                          % space for reverse index
jfint(ifint)=1:mn;                          % create reverse index
iy=1:m;                                     % index of output values
kfint=jfint(iy)-iy;                         % next lower x knot for each y knot (in range 0:n)
nle=jfint(m+1)-1;                           % number of kfint=0 values at start (these need extrapolation)
nhe=mn-jfint(mn);                           % number of kfint=n values at end (these need extrapolation)
% could check for errors: nle+nhe=m means no interpolation possible
cl=1+abs(kfint-1);                          % lower x knot for interpolation
ch=n-abs(kfint-n+1);                        % upper x knot for interpolation
frac=(y-x(cl))./(x(ch)-x(cl));              % fraction of the way that y knot is from x(cl) to x(ch)
vv=[1-frac(:) frac(:)];                     % interpolation weights (assuming extrapolation at each end)
switch ws                                   % check extrapolation method at low end
    case 0                                  % set interpolation weights to zero for first nle y knots
        vv(1:nle,:)=0;
    case 1                                  % set first nle values of y to x(1)
        vv(1:nle,:)=repmat([0 1],nle,1);
end
switch we                                   % check extrapolation method at high end
    case 0                                  % set interpolation weights to zero for  last nhe y knots
        vv(m-nhe+1:m,:)=0;
    case 1                                  % set last nhe values of y to x(1)
        vv(m-nhe+1:m,:)=repmat([1 0],nhe,1);
end
inm=sparse(repmat(iy,1,2),[cl ch],vv(:),m,n);   % create m x n interpolation matrix
if ~nargout
    for i=1:m
        si=sum(inm(i,:),2);
        msk=inm(i,:)~=0;
        sm=sum(msk);                        % number of non-zero weights for y(i)
        if sm>0
            plot(y(i),si,'ob',x(msk),inm(i,msk),'+b',[x(msk); repmat(y(i),1,sm)],[inm(i,msk); repmat(si,1,sm)],'-b');
        else
            plot(y(i),si,'ok')
        end
        hold on
    end
    axisenlarge(-1.05);
    xlim=get(gca,'xlim');
    plot(xlim,[0 1; 0 1],':k');
    hold off
end