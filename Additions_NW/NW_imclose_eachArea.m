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

MARGIN=10;

outputIm=inputIm;

%% find indices of areas (typically 0,1,2,...)
AreaNrs=unique(inputIm)';  % need row vector
idx=find(AreaNrs~=0); % background=0
AreaNrs=AreaNrs(idx);
clear idx


% find integers (via modulo because 'unique' has double as output???)
% idxnotint=find(mod(AreaNrs,1))~=0;
%if ~isempty(idxnotint)
if any(mod(AreaNrs,1))
    error('Some areas have non-integer values!')
end

%% loop over all areas and fill fuzzy edges (imclose)
%tic
sizeInputIm=size(inputIm); % MW
for runidx=AreaNrs   
   
    %%
    % The 3 commented lines below are shorter in script, but take 2X as
    % long as code below
    %{
    SubAreaonlyIm=(inputIm==runidx);
    SubAreaonlyIm=imclose(SubAreaonlyIm,struct_el);
    outputIm(SubAreaonlyIm==1)=runidx;
    %}
    
     
    
    %%    
    % make a little image containig the cell of interest
    [i,j]=ind2sub(sizeInputIm,find(inputIm==runidx));
    toSelecti = [max(1,min(i)-MARGIN):min(max(i)+MARGIN,sizeInputIm(1))];
    toSelectj = [max(1,min(j)-MARGIN):min(max(j)+MARGIN,sizeInputIm(2))];    
    SubAreaonlyIm = inputIm(toSelecti, toSelectj);
    SubAreaonlyIm = (SubAreaonlyIm==runidx);
    
    % close it
    SubAreaonlyIm=imclose(SubAreaonlyIm,struct_el);
    
    % update output image accordingly
    for idxi = 1:numel(toSelecti)
        for idxj = 1:numel(toSelectj)
            
            i=toSelecti(idxi);
            j=toSelectj(idxj);
            
            if SubAreaonlyIm(idxi,idxj)
                outputIm(i, j) = runidx;
            end
            
        end
    end
    

end
%toc

% figure(); imshow(outputIm,[]);

