function [A, Nn, Nrg, Ff]=gciFourierAnalysis(sp,fs,gci)

wname=@rectwin;
gci(end+1)=length(sp);

Nf=max(diff(gci));
Nt=length(gci)-1;
A=nan(Nt,Nf);
Ff=nan(Nt,Nf);
Nn=zeros(1,Nt);
Nrg=ones(1,Nt);

for ik=1:Nt
    spseg=sp(gci(ik):gci(ik+1)-1);
    %Nrg(ik)=sqrt(spseg(:)'*spseg(:));
    spseg=spseg/Nrg(ik);
    N=length(spseg);
    Nn(ik)=N;
    A(ik,1:N)=fft(spseg);
    Ff(ik,1:N)=(fs/N)*(0:(N-1));
    
    
    %Nfft=N;
    %A(ik,1:Nfft)=fft(spseg.*window(wname,N),Nfft);
    %Ff(ik,1:Nfft)=(fs/Nfft)*(0:(Nfft-1));
   
end

%A(ik,1:Nfft)=2*pi*fft(spseg.*window(wname,N),Nfft)/N;