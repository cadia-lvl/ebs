function out=mymap(fun, par, list1, list2)

%  Run a function on lists of files
%
%         out=mymap(fun, par, list1, list2)
%
%   Inputs:     
%       fun         A function handle that is run as out=feval(@fun,file1,file2,par)
%       par         Structure as an input for @fun
%       list1       List of files used for input to @fun
%       list2       (Optional) List of files used for output for @fun
%
%   Outputs written in outfile:
%       out         A structure with outputs from @fun, if provided
%
%   References:  http://en.wikipedia.org/wiki/Map_(higher-order_function)
% similar to cellfun
%
%    NB.  Could make more general for n lists using varargin etc.
%
%
%**************************************************************************
% Author:           J. Gudnason
% Date:             16 March 2009
% $Revision: 1.3 $
% $Date: 2009/09/17 08:19:57 $
%**************************************************************************


no=nargout(fun);
cnt1=countLines(list1);
if nargin>3
    cnt2=countLines(list2);
    if cnt1~=cnt2
        error('The filelists have different number of lines');
    end;
    fid2=fopen(list2,'r');
end;
fid1=fopen(list1,'r');

out=cell(cnt1,no);
h = waitbar(0,['Running  ' func2str(fun) ' on ' list1]);
for ir=1:cnt1
    file1=fgetl(fid1);
    disp(file1);
    try
        if nargin>3
            file2=fgetl(fid2);
            [a{1:no}] = feval(fun,file1,file2,par);
        else
            [a{1:no}] = feval(fun,file1,par);
        end;
    catch
        %error(['In processing file ' file1 '.']);
        [a{1:no}]=file1;%warning(file1);
    end
    if no
        out(ir,:)=a;
    end;
    waitbar(ir/cnt1,h);
end;
close(h);
fclose(fid1);
if nargin>3, fclose(fid2); end;


