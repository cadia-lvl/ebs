function timeMarkers = vuvExtract(ifile, praatScriptDir, praatCmd)

%VUVEXTRACT is a wrapper function for the Praat script vuvExtract.praat [timeMarkers] = (ifile, praatScriptDir, praatCmd)
%
%  Inputs:  ifile           is the input wav file
%           praatScriptDir  is the directory of the vuvExtract.praat
%           praatCmd        is the Praat command to be executed with the with the 'system' command
%
% Outputs:  timeMarkers(Nx2)    matrix defining the voiced segments (first
%                               column is the start and second the end of a 
%                               voiced segment (in seconds))


if nargin < 3
    comp = computer;
    
    if strcmp(comp, 'MACI64')
        praatCmd = '/Applications/Praat.app/Contents/MacOS/Praat';
    elseif strcmp(comp, 'PCWIN64')
        praatCmd = 'C:\Users\Lalli\Documents\Rannsoknir\praatcon.exe';
    else
        error('Dont know where Praat is.');
    end;
end

vuvPraatScript = [praatScriptDir '/vuvExtract.praat'];
%vuvPraatScript='/Users/jg/Work/PraatTools/vuvExtract.praat';
%vuvPraatScript='/Users/jg/GoogleDrive/Work/Software/Praat/vuvExtract.praat';
timeStep = 0.01;
pitchFloor = 60.0;
maxNoCand = 15;
veryAccurate = 'no';
silenceTh = 0.03;
voicingTh = 0.7;
octiveCost = 0.01;
octiveJumpCost = 0.35;
vuCost = 0.14;
pitchCeiling = 400.0;

cmd= [praatCmd ' ' vuvPraatScript ' ' ifile ' tmp.txt ' ...
    num2str(timeStep) ' ' ...
    num2str(pitchFloor) ' ' ...
    num2str(maxNoCand) ' ' ...
    veryAccurate ' ' ...
    num2str(silenceTh) ' ' ...
    num2str(voicingTh) ' ' ...
    num2str(octiveCost) ' ' ...
    num2str(octiveJumpCost) ' ' ...
    num2str(vuCost) ' ' ...
    num2str(pitchCeiling)];

[status, result] = system(cmd);

if ~isempty(result)
    warning(result);
end

fid = fopen([praatScriptDir '/tmp.txt'],'r');
ii=1;
st2 = '0';
st1 = '0';
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if strcmp(tline,'"V"')
        timeMarkers(ii,1) = str2num(st2);
        timeMarkers(ii,2) = str2num(st1);
        ii = ii+1;
    end
    st2 = st1;
    st1 = tline;
end
fclose(fid);
%!cat /Users/jg/GoogleDrive/Work/Software/Praat/tmp.txt
delete([praatScriptDir '/tmp.txt']);

