function A = changem_clone(A, newCode, oldCode)
%function A = changem_clone(A, newCode, oldCode)
% Function that mimics matlab's changem() function, but does not require
% the MAP toolbox.
%
% Replaces instances of oldCode(i) in matrix A with newCode(i).
% e.g.
% A = changem_clone([1,1,2,3;1,1,2,3], [7,8], [1,2])
% % Result: A = [7 7 8 3;7 7 8 3]
% 
% Thanks to Stackoverflow user Rody Oldenhuis
% http://stackoverflow.com/questions/13812656/elegant-vectorized-version-of-changem-substitute-values-matlab

[a,b] = ismember(A,oldCode);

A(a) = newCode(b(a));