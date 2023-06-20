function SNR=calculateSNR(s,vsr)

nsu=length(vsr);                                
ssq=sum(s(1:nsu).^2);
SNR=-db(sum((vsr(1:nsu)-s(1:nsu)).^2)/ssq);