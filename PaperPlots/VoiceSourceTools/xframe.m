function [X,dT,nrg]=xframe(sp,fs,frameStart,Nsamp)


%XFRAME Create data matrix from a signal and frame boundaries [X,dT,nrg]=(sp,fs,frameStart,Nsamp)
%
%  Inputs:  sp              The input signal
%           fs              The sampling frequency
%           frameStart      The M indices of each starting sample of the frames in the signal
%           Nsamp           The number of samples each frame is resampled to
%
% Outputs:  X               The data matrix (MxNsamp)
%           dT              The time length of each frame
%           nrg             The energy of each window (sqrt(x.'*x))


numFrames=length(frameStart);
X=zeros(numFrames,Nsamp);
dT=zeros(1,numFrames);
nrg=zeros(1,numFrames);
frameStart(end+1)=length(sp);
for ii=1:numFrames
    NSeg=frameStart(ii+1)-frameStart(ii)+1;
    spSeg=sp(frameStart(ii):frameStart(ii+1));
    [x]=resample(spSeg,Nsamp,NSeg);
    nrg(ii)=sqrt(x.'*x);
    X(ii,:)=x./(nrg(ii)+eps);
    dT(ii)= NSeg/fs;
end