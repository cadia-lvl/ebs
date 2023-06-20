addpath ../utils/
addpath ../tools/

par=projParam();

sstr=['test/dr3/mjmp0/sx'];
disp(sstr);
makeListsTIMIT(par.pth,sstr);
out{1,1}=mymap(@stftMelCopySynth,par.db,par.pth.splistname);


