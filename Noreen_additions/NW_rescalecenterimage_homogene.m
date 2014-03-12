function imagescaled = NW_rescalecenterimage_homogene(origimage,rescalefactor)
% rescales image by factor 'rescalefactor'. Then centers it so that new
% central point is the same as the original center. Then image is cut to
% have original dimensions

imagedim=size(origimage);   % size of rescaled origimage must be cut to this initial size
              
% find centre points of image (not yet rounded) and their
% shift after rescaling
centerpoints=0.5*size(origimage);
newcenterpoints=centerpoints*rescalefactor;
centerpointsshift=round(newcenterpoints-centerpoints);
rowshift=centerpointsshift(1); columnshift=centerpointsshift(2); % only for better readability;
              
%rescale origimage
imageintermed=imresize(origimage,rescalefactor,'bicubic'); % bicubic should be default anyway.
                    % leave imresize. probably also fine with imresize_old (not tested)
                    % I think bicubic is best interpolaiton ... (NW 2012)
% fit in old frame
imagescaled(1:imagedim(1),1:imagedim(2))=imageintermed(1+rowshift:imagedim(1)+rowshift,1+columnshift:imagedim(2)+columnshift);
                        % how to choose the shift: imaginge scale factor=1, then rowshift, columnshift=0 
                        % and first coordinate must be mapped on first coordinate

                        % for some reason, this leads to the fact that (when scaling by a
                        % factor 2) the new center is a 2x2 mattrix with oldcentercoordinates and
                        % -1 in each coordinate, e.g. (45,50) -> (44 45, 49 50)
                        % should not be a big problem because origimage
                        % still will be shifted

