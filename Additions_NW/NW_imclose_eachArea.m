function outputIm=NW_imclose_eachArea(inputIm,struct_el)
% This function performs the morphological operation 'imclose' on each
% subarea (=region with a specific integer number) of inputIm.
% Function is called from NW_segphase_richMed
%
% INPUT:
% inputIm: image in which each area has a different integer value
%          (background=0). Typical input: segmentation image or in-between steps
%          during segmentation.
%          Note: If input image is binary, it can be converted to different
%          integer per connected region via ImInt=bwlabel(ImBin)
% struc_el: structuring element for morphological operation. 
%           E.g. strel('disk',5)
%

outputIm=inputIm;

% find indices of areas (typically 0,1,2,...)
AreaNrs=unique(inputIm)';  % need row vector
idx=find(AreaNrs~=0); % background=0
AreaNrs=AreaNrs(idx);
clear idx
% find integers (via modulo because 'unique' has double as output???)
idxnotint=find(mod(AreaNrs,1))~=0;
if ~isempty(idxnotint)
    error('Some areas have non-integer values!')
end

%loop over all areas and fill fuzzy edges (imclose)
for runidx=AreaNrs
    SubAreaonlyIm=(inputIm==runidx);
    SubAreaonlyIm=imclose(SubAreaonlyIm,struct_el);
    outputIm(SubAreaonlyIm==1)=runidx;
end
