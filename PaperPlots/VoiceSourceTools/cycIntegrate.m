function u = cycIntegrate(uu, fs, gci)

tsrat=0.6;
terat=1-tsrat;
T=gci(2)-gci(1);

u=zeros(size(uu));
nn=1:(gci(1)-ceil(tsrat*T));
uuseg=uu(nn);
useg = cumsum(uuseg)/fs;
useg=useg-linspace(0,useg(end),length(useg))';  % Enforce equal pressure at end
u(nn)=useg;

for ig=1:length(gci)
    nn=(-ceil(tsrat*T):ceil(terat*T))+gci(ig);
    % If the indices exceed the signal boundaries then cut (and adjust T)
    nn=nn(nn>0); 
    nn=nn(nn<length(u));
    T=length(nn)-2;    % T should be 2 samples short of length nn

    
    uuseg=uu(nn);
    useg = cumsum(uuseg)/fs;
    useg=useg-linspace(0,useg(end),length(useg))';  % Enforce equal pressure at end
    u(nn)=useg;
    
    T=gci(min(ig+1,length(gci)))-gci(ig);
end;

nn=(gci(end)+ceil(terat*T)):length(u);
uuseg=uu(nn);
useg = cumsum(uuseg)/fs;
useg=useg-linspace(0,useg(end),length(useg))';  % Enforce equal pressure at end
u(nn)=useg;