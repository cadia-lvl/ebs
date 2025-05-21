function name = computername()
% COMPUTERNAME returns the name of the computer (hostname)
% name = getComputerName()

% See also SYSTEM, GETENV, ISPC, ISUNIX
%
% m j m a r i n j (AT) y a h o o (DOT) e s
% (c) MJMJ/2007
%
persistent nm
if ~numel(nm)
[ret, name] = system('hostname');   
if ret ~= 0,
   if ispc
      name = getenv('COMPUTERNAME');
   else      
      name = getenv('HOSTNAME');      
   end
end
nm=strtrim(name);  % remove any junk at the end
end
name=nm;


