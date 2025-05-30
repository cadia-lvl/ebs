function PatchUnevenSpectrogram(fh,magA,fs,gci,Nn)

for ik=1:length(Nn)
    tb=gci(ik)/fs;
    te=gci(ik+1)/fs;
    fb=0; 
    fe=fs/Nn(ik);
    for is=1:ceil(Nn(ik)/2)
        patch(1000*[tb te te tb],[fb fb fe fe]/1000,magA(ik,is),'EdgeColor','none')
        fb=fb+fs/Nn(ik);
        fe=min(fe+fs/Nn(ik),fs/2);
    end
end
cb=colorbar;
cb.Label.String='Power in dB';
xlabel('Time [ms]')
ylabel('Frequency [kHz]')