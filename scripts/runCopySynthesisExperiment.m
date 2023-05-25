addpath ../utils/
addpath ../tools/

par=projParam();

traintest={'train','test'};
for it=1%:2   %Uncomment to run on test as well
    for idr=1%:8 %Uncomment to run DR2-8 as well
        sstr=[traintest{it} '/dr' num2str(idr)];
        disp(sstr);
        makeListsTIMIT(par.pth,sstr);
        out{it,idr}=mymap(@stftMelCopySynth,par.db,par.pth.splistname);
        save('RawRes.mat','out','it','idr','-v7.3');
    end
end;
