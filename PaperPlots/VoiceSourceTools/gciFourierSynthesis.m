function sp=gciFourierSynthesis(A,Nn,Nrg)

wname=@rectwin;

%Nf=max(diff(gci));
Nt=size(A,1);

sp=zeros(sum(Nn),1);
n=1;

for ik=1:Nt
    sp(n:(n+Nn(ik)-1))=Nrg(ik)*ifft(A(ik,1:Nn(ik)),'symmetric');
    n=n+Nn(ik);
end