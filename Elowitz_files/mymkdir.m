function [status,msg,id] = mymkdir(absoluteDir);
% This is a wrapper around mkdir, because mkdir with one argument 
% behaves differently in matlab 6.5 vs. matlab 7.0.  This behaves
% like matlab 7.0 mkdir when called with one argument, and it assumes 
% the one argument is an absolute directory.  Works on Win & Unix.

v = version;
if (v(1)=='7')
  [status, msg, id] = mkdir(absoluteDir);
  return
end

% below, we assume we're in version matlab 6.*

if isunix
  [status, msg, id] = mkdir('/',absoluteDir);
  return
end

if ispc
  [status, msg, id] = mkdir(absoluteDir(1:2),absoluteDir(3:end));
  return
end

error('mymkdir cant figure out what OS you have!')

