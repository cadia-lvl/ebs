function stftr=mel2spec_nzp(stftm,fs,nfft,framelen,preserveDC,stftangle)


[numMelFilters,nframes]=size(stftm);


stftr=zeros(nfft,nframes);

for ii=1:nframes
    nfftp=1+floor(framelen(ii)/2);  
    %nfftn=framelen(ii)-nfftp; % Number of negative frequency cofficients
    if preserveDC
        mbm=v_filtbankm(numMelFilters-1,framelen(ii),fs,0,fs/2,'m');
        mbm=vertcat(sparse(1,1,1,1,nfftp),mbm); % preappend an extra row to preserve the DC value
    else
        mbm=v_filtbankm(numMelFilters,framelen(ii),fs,0,fs/2,'m');
    end
    imbm=pinv(full(mbm));
    posFr=sqrt(max(imbm*stftm(:,ii),0)).*exp(-1i*stftangle(1:nfftp,ii));  %Includes DC and Nyquist for even lengths

    if rem(framelen(ii),2)==0  % Even length sequence
        negFr=conj(posFr((nfftp-1):-1:2));  % The Nyquist bin is included in posFr so ommitted here
    else  % Odd length sequence
        negFr=conj(posFr(nfftp:-1:2));    % There is no Nyquist bin for odd sequences
    end;

    stftr(1:framelen(ii),ii)=[posFr; negFr]';
end
