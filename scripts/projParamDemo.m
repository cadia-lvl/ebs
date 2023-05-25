function par=projParamDemo()

par.opt.listopt='sfp'; %(s = make speech file list, x = make lx list, f = make file output list, f = make phn list)
par.opt.fetext='mat';
par.opt.format='timit'; %('vpd', 'nist', 'wav' etc.  not all supported yet.)

if strcmp(par.opt.format,'timit')
    fs=16000;
else
    fs=20000;
end

%Parameters needed for sdftMelCopySynthDB
par.db.frameMethod = {'Fixed', 'Epoch', 'EpochAdj'};
par.db.gcifrac=0.3;                                        % position of GCI in analysis frame
par.db.maxoff=10;                                          % maximum offset in samples (obsolete)
par.db.maxtoff=0.625e-3;                                   % maximum offset in seconds
par.db.numMels=[2 5 10 15 20 30 40 60 80 100 140];
par.db.pitchlim=[40 50 400];                               % {min targt max} pitch (Hz)
par.db.periodlimk=round(fs./par.db.pitchlim);              % {max targt mmin} pitch periods (samples)
par.db.nfft=2*ceil(par.db.periodlimk(1)/2);                % length of fft (force to be an even number)
par.db.fs=fs;
par.db.doPESQ=1;
par.db.PESQpath="<local computer full path for PESQ>/tools/PESQ";
par.db.GCImethod='YAGA';                                   % Can be 'YAGA' or 'SEDREAMS'
par.db.MELmethod='filtbankm';                              % Can be 'melbankm' or 'filtbankm'

%Parameters needed for making filelists /globbing
par.pth.speechpth = '<local computer full path for TIMIT>';
par.pth.frm='wav';
par.pth.fetpth = ['./Features'];
par.pth.splistname = 'splist.scp';
par.pth.lxlistname = 'lxlist.scp';
par.pth.fetlistname = 'fetlist.scp';
par.pth.phlistname = 'phlist.scp';
par.pth.fileinfo = 'fileinfo.mat';
