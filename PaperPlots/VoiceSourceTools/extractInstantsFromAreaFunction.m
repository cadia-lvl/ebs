function [gci, goi] = extractInstantsFromAreaFunction(AreaFnc,th)

if nargin<2
    %AreaFnc is in cm^2
    %thR=0.01;
    
    %AreaFnc is in m^2
    thR=1e-6;
else
    thR=th;
end

lA=length(AreaFnc);
sA=round(0.2*lA);
eA=round(0.8*lA);

AreaFnc=AreaFnc-min(AreaFnc(sA:eA));

inst=diff(AreaFnc>thR);

gci=find(inst==-1);
goi=find(inst==+1);

%Assumtion that gci and goi occur in tandem, just force gci to be first and
%goi to be last in sequence:

goi(goi < gci(1))=[];
gci(gci > goi(end))=[];

% Remove ones close to beginning
ix=median(diff(gci))<gci;
iy=(length(AreaFnc)-median(diff(gci)))>goi;
gci=gci(and(ix,iy))';
goi=goi(and(ix,iy))';