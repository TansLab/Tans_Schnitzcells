function PN_segLoop(p, array_segRange, array_array_edgeSlices, array_array_fillingEdgeSlices, array_minNumEdgePixels, array_botHatSize, removeSmallCells) 

%
if ~exist('removeSmallCells')
  removeSmallCells = 1;
end

% ANALYSIS
for i0 = array_botHatSize
  for i1 = array_minNumEdgePixels
    for i2 = array_segRange
      for i3 = [1:size(array_array_edgeSlices,1)]
        for i4 = [1:size(array_array_fillingEdgeSlices,1)]
          i3Nrs = array_array_edgeSlices(i3,:);
          i3Nrs(i3Nrs==0) = [];
          i3String = num2str(i3Nrs);
          i3String(i3String==' ') = [];

          i4Nrs = array_array_fillingEdgeSlices(i4,:);
          i4Nrs(i4Nrs==0) = [];
          i4String = num2str(i4Nrs);
          i4String(i4String==' ') = [];


          directory = ['improved_edge' i3String '_fill' i4String '_pix' str3(i1) '_hat' num2str(i0) '\seg' str3(i2) extra_string_filename '\']; 

          PN_segmoviephase(p, 'minCellLengthConservative',20,'removeSmallCells',removeSmallCells,'segRange',i2,'edgeSlices',i3Nrs,'minNumEdgePixels',i1,'fillingEdgeSlices',i4Nrs,'botHatSize',i0,'DJK_saveDir',directory); 

        end
      end
    end
  end
end
