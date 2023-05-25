function phn=extractPhn(phnfile, fetfile)

%   EXTRACTPHN extract phone information from TIMIT datafiles and add to mat file
%
%   extractPhn(phnfile, fetfile, ph)
%
%   Inputs:
%       phnfile    TIMIT phone file
%       fetfile    Mat file for outputs (field phn with data added)
%
%   References:
%
%**************************************************************************
% Author:           J. Gudnason 
% Date:             May 2009
% Changed :         August 2020
%********************************************************************

if nargin<3
    ph=timitPhones();
end
fid = fopen(phnfile, 'r');
txt = textscan(fid, '%f %f %s\n');
fclose(fid);
phn.tstart = txt{1};
phn.tend = txt{2};
phn.name = txt{3};

for ii=1:length(phn.name)
    phn.in(ii)=find(strcmpi(phn.name{ii},ph));
end

phn.vus=timitVUS(phn);

if nargout==0
    if exist(fetfile,'file')
        save(fetfile,'phn','-append');
    else
        save(fetfile,'phn');
    end
end

function vus=timitVUS(phn)

vIn=[19:25 26:32 33:52];
sIn=[53:57];

n=phn.tend(end);
segN=length(phn.tend);
vus=2*ones(1,n);

for is=1:segN
    if any(phn.in(is)==sIn)
        vus(phn.tstart(is)+1:phn.tend(is))=1;
    end
    if any(phn.in(is)==vIn)
        vus(phn.tstart(is)+1:phn.tend(is))=3;
    end
end


function ph=timitPhones()

%Stops:   
ph{1}='b';          %bee           BCL B iy
ph{2}='d';          %day           DCL D ey
ph{3}='g';          %gay           GCL G ey 
ph{4}='p';          %pea           PCL P iy
ph{5}='t';          %tea           TCL T iy
ph{6}='k';          %key           KCL K iy
ph{7}='dx';         %muddy, dirty  m ah DX iy, dcl d er DX iy
ph{8}='q';          %bat           bcl b ae Q
%Affricates:
ph{9}='jh';         %joke          DCL JH ow kcl k
ph{10}='ch';         %choke         TCL CH ow kcl k
%Fricatives:
ph{11}='s';          %sea           S iy
ph{12}='sh';         %she           SH iy
ph{13}='z';          %zone          Z ow n
ph{14}='zh';         %azure         ae ZH er
ph{15}='f';          %fin           F ih n
ph{16}='th';         %thin          TH ih n
ph{17}='v';          %van           V ae n 
ph{18}='dh';         %then          DH e n
%Nasals:
ph{19}='m';          %mom           M aa M
ph{20}='n';          %noon          N uw N
ph{21}='ng';         %sing          s ih NG
ph{22}='em';         %bottom        b aa tcl t EM
ph{23}='en';         %button        b ah q EN
ph{24}='eng';        %washington    w aa sh ENG tcl t ax n
ph{25}='nx';         %winner        w ih NX axr
%Semivowels and Glides:
ph{26}='l';          %lay           L ey
ph{27}='r';          %ray           R ey
ph{28}='w';          %way           W ey
ph{29}='y';          %yacht         Y aa tcl t
ph{30}='hh';         %hay           HH ey
ph{31}='hv';         %ahead         ax HV eh dcl d
ph{32}='el';         %bottle        bcl b aa tcl t EL
%Vowels:
ph{33}='iy';         %beet          bcl b IY tcl t
ph{34}='ih';         %bit           bcl b IH tcl t 
ph{35}='eh';         %bet           bcl b EH tcl t
ph{36}='ey';         %bait          bcl b EY tcl t
ph{37}='ae';         %bat           bcl b AE tcl t
ph{38}='aa';         %bott          bcl b AA tcl t
ph{39}='aw';         %bout          bcl b AW tcl t
ph{40}='ay';         %bite          bcl b AY tcl t
ph{41}='ah';         %but           bcl b AH tcl t
ph{42}='ao';         %bought        bcl b AO tcl t
ph{43}='oy';         %boy           bcl b OY
ph{44}='ow';         %boat          bcl b OW tcl t
ph{45}='uh';         %book          bcl b UH kcl k
ph{46}='uw';         %boot          bcl b UW tcl t
ph{47}='ux';         %toot          tcl t UX tcl t
ph{48}='er';         %bird          bcl b ER dcl d
ph{49}='ax';         %about         AX bcl b aw tcl t
ph{50}='ix';         %debit         dcl d eh bcl b IX tcl t
ph{51}='axr';        %butter        bcl b ah dx AXR
ph{52}='ax-h';       %suspect       s AX-H s pcl p eh kcl k tcl t
%Others:
ph{53}='pau';     %pause
ph{54}='epi';     %epenthetic silence
ph{55}='h#';      %begin/end marker (non-speech events)
ph{56}='1';       %primary stress marker
ph{57}='2';       %secondary stress marker
%Closure intervals of stops
ph{58}='bcl';
ph{59}='dcl';
ph{60}='gcl';
ph{61}='pcl';
ph{62}='tck';
ph{63}='kcl';
ph{64}='tcl';

