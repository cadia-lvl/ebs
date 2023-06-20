function sfr=gs_frameExtLPC(s,framex,par)

Nmax=par.nfft;
p=16;

nframe=length(framex); % number of frames
framekk=[[1 framex(1:end-1)+1];framex];                         %  start and end samples for each frame
framelen=framekk(2,:)-framekk(1,:)+1;
if max(framelen)>Nmax
    error('Maximum framelength is larger than Nmax');
elseif min(framelen)<2*p
    error('Minimum frame length is smaller than 2*p');
end;

sfr=zeros(Nmax,nframe);

for i=1:nframe
    sseg=s(framekk(1,i):framekk(2,i));
    sfr(1:framelen(i),i)=sseg;
    if framelen(i)<Nmax
        
        ar=v_lpcauto(sseg,p,[],'r','e');
        uu=filter(ar,1,sseg);
        [~,zf]=filter(1,ar,uu);  %First output should be = sseg
        sExt=filter(1,ar,zeros((Nmax-framelen(i)),1),zf);

        sfr(framelen(i)+1:Nmax,i)=sExt;

    end;
end;
