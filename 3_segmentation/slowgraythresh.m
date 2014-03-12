function level = slowgraythresh(I)
% function level = slowgraythresh(I)
%
% matlab routine changed to handle uint16's
% 
%GRAYTHRESH Compute global image threshold using Otsu's method.
%   LEVEL = GRAYTHRESH(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0, 1].
%   GRAYTHRESH uses Otsu's method, which chooses the threshold to minimize
%   the intraclass variance of the thresholded black and white pixels.
%
%   Class Support
%   -------------
%   The input image I can be of class uint8, uint16, or double.  LEVEL
%   is a double scalar.
%
%   Example
%   -------
%       I = imread('blood1.tif');
%       level = graythresh(I);
%       BW = im2bw(I,level);
%       imshow(BW)
%
%   See also IM2BW.

%   Copyright 1993-2001 The MathWorks, Inc.  
%   $Revision: 1.1.1.1 $  $Date: 2004/12/08 19:10:32 $

% Reference:
% N. Otsu, "A Threshold Selection Method from Gray-Level Histograms,"
% IEEE Transactions on Systems, Man, and Cybernetics, vol. 9, no. 1,
% pp. 62-66, 1979.

% One input argument required.
error(nargchk(1,1,nargin));

if ~isa(I,'uint16'),
    disp('sorry -- this special version was created only for uint16''s -- tough luck, man...');
    return;
end;

if (~isa(I,'uint8') & ~isa(I,'double') & ~isa(I,'uint16'))
    error('Input image must be uint8, uint16, or double.')
end

% Convert all N-D arrays into a single column.  Convert to uint8 for
% fastest histogram computation.



% I = im2uint8(I(:));

num_bins = 65536;
counts = imhist(I,num_bins);

% Variables names are chosen to be similar to the formulas in
% the Otsu paper.
p = counts / sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:num_bins)');
mu_t = mu(end);

% Save the warning state and disable warnings to prevent divide-by-zero
% warnings.
state = warning;
warning off;
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

% Restore the warning state.
warning(state);

% Find the location of the maximum value of sigma_b_squared.
% The maximum may extend over several bins, so average together the
% locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
% then return 0.
maxval = max(sigma_b_squared);
if isfinite(maxval)
    idx = mean(find(sigma_b_squared == maxval));
    
    % Normalize the threshold to the range [0, 1].
    level = (idx - 1) / (num_bins - 1);
else
    level = 0.0;
end
