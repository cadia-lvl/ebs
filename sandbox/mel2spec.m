function stftr=mel2spec(stftm,fs,nfft,preserveDC,stftangle)


numMelFilters=size(stftm,1);
nfftp=1+floor(nfft/2);
if preserveDC
    mbm=v_filtbankm(numMelFilters-1,nfft,fs,0,fs/2,'m');
    mbm=vertcat(sparse(1,1,1,1,nfftp),mbm);
else
    mbm=v_filtbankm(numMelFilters,nfft,fs,0,fs/2,'m');
end;
imbm=pinv(full(mbm));

stftr=sqrt(max(imbm*stftm,0)).*exp(1i*stftangle(1:nfftp,:));
stftr=[stftr; conj(stftr(end-1:-1:2,:))];
