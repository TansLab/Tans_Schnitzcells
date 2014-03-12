function level = joeslowgraythresh(I,name)
% function level = joeslowgraythresh(I,name)
%
% Joe's revised version of slowgraythresh... 
%   History: 
%     - slowgraythresh was matlab graythresh changed to handle uint16's
%     - joe changed to simpler approach, set threshold to 1 sigma below mean 
%   optional input argument "name" turns on figure output including title "name"
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
%   $Revision: 1.2 $  $Date: 2005/01/12 02:09:10 $

% Reference:
% N. Otsu, "A Threshold Selection Method from Gray-Level Histograms,"
% IEEE Transactions on Systems, Man, and Cybernetics, vol. 9, no. 1,
% pp. 62-66, 1979.

% One input argument required.
error(nargchk(1,2,nargin));

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

% JCR: debugging threshold problems... 
    
if (nargin == 2)

  figure
  subplot(3,1,1)
  plot([1:10000]/num_bins,counts(1:10000),'b-','MarkerSize',2);
  hold on
  xlabel('pixel gray value');
  ylabel('pixel counts');
  h = line([level level],[0 max(counts(1:10000))]);
  set(h,'Color',[0 0 1]);
  title(['Threshold computed for ' name]);
  
  subplot(3,1,2)
  plot([1:10000]/num_bins,sigma_b_squared(1:10000),'r-')
  hold on
  xlabel('pixel gray value');
  ylabel('\sigma_b ^2');
  h = line([level level],[0 maxval]);
  set(h,'Color',[0 0 1]);
  
  subplot(3,1,3)
  d_sigma_b_squared = diff(sigma_b_squared(1:10000));
  max_d = max(d_sigma_b_squared);
  plot([2:10000]/num_bins,d_sigma_b_squared,'r-')
  hold on
  grid on
  xlabel('pixel gray value');
  ylabel('1st derivative of \sigma_b ^2');
  h = line([level level],[0 max_d]);
  set(h,'Color',[0 0 1]);

end

% JCR: For bacillus, this threshold seems to be working better across frames
pixels = double(reshape(I,prod(size(I)),1))/65535;
[muhat,sigmahat] = STnormfit(pixels);  %%%change SJT
% hack- replace muhat with median!
muhat = median(pixels);
% another hack- for frames having mostly background, default level seems ok
if muhat < level
  % Try only fixing the level for when the default level is broken
  level = muhat-sigmahat*0.5;
end

if (nargin == 2)
  
  subplot(3,1,1)
  h = line([muhat muhat],[0 max(counts(1:10000))]);
  set(h,'Color',[0 1 0]);
  h = line([muhat-sigmahat muhat-sigmahat],[0 max(counts(1:10000))]);
  set(h,'Color',[0 1 0]);
  h = line([muhat+sigmahat muhat+sigmahat],[0 max(counts(1:10000))]);
  set(h,'Color',[0 1 0]);
  h = line([level level],[0 max(counts(1:10000))]);
  set(h,'Color',[1 0 1]);

  subplot(3,1,2)
  h = line([level level],[0 maxval]);
  set(h,'Color',[1 0 1]);

  subplot(3,1,3)
  h = line([level level],[0 max_d]);
  set(h,'Color',[1 0 1]);

end
