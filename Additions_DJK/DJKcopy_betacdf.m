function p = DJKcopy_betacdf(x,a,b);
%BETACDF Beta cumulative distribution function.
%   P = BETACDF(X,A,B) returns the beta cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINC does the computational work.
%
%   See also BETAFIT, BETAINV, BETALIKE, BETAPDF, BETARND, BETASTAT, CDF,
%            BETAINC.
   
%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.5.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.9.2.6 $  $Date: 2004/12/24 20:46:45 $

if nargin<3, 
   error('stats:betacdf:TooFewInputs','Requires three input arguments.'); 
end

[errorcode x a b] = DJKcopy_distchck(3,x,a,b);

if errorcode > 0
   error('stats:betacdf:InputSizeMismatch',...
         'Requires non-scalar arguments to match in size.');
end

% Initialize P to 0.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

p(a<=0 | b<=0) = NaN;

% If is X >= 1 the cdf of X is 1. 
p(x >= 1) = 1;

k = find(x > 0 & x < 1 & a > 0 & b > 0);
if any(k)
   p(k) = betainc(x(k),a(k),b(k));
end

% Make sure that round-off errors never make P greater than 1.
p(p > 1) = 1;
