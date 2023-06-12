function par=projParamDemo()

par.pth.listopt='s'; %(s = make speech file list, x = make lx list, f = make file output list, f = make phn list)
par.pth.fetext='mat';
fs=16000;

%Parameters needed for making filelists /globbing
par.pth.speechpth = '/Users/username/TIMIT/';
par.pth.frm='wav';
par.pth.fetpth = '';
par.pth.splistname = 'splist.scp';
par.pth.fetlistname = 'fetlist.scp';
par.pth.phlistname = 'phlist.scp';

%Parameters needed for sdftMelCopySynthDB
par.db.frameMethod = {'Fixed', 'Epoch', 'EpochAdj'};
par.db.gcifrac=0.3;                                        % position of GCI in analysis frame
par.db.maxoff=10;                                          % maximum offset in samples
par.db.numMels=[2 5 10 15 20 30 40 60 80 100 140];
par.db.fbankmethod='filtbankm';
par.db.preserveDC=1;
par.db.pitchlim=[40 50 400];                               % {min targt max} pitch (Hz)
par.db.periodlimk=round(fs./par.db.pitchlim);                     % {max targt mmin} pitch periods (samples)
par.db.nfft=2*ceil(par.db.periodlimk(1)/2); % length of fft (force to be an even number)
par.db.fs=fs;
par.db.doPESQ=1;
par.db.PESQpath='../utils/PESQ';
par.db.GCImethod='YAGA';