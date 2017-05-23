% Also stolen from NW_tracker_centroid_vs_area

function [Lc_fullsize_centered, Lc_fullsize, Lc] = MW_loadLcData(segFile);
    load(segFile);
    Lc_fullsize           = zeros(phaseFullSize);
    Lc_fullsize_centered  = zeros(phaseFullSize);

    % if segmentation is not approved, get unapproved segmentation 
    if ~(exist('Lc') == 1) 
      Lc = LNsub;
    end

    % put Lc back in original location in fullsize
    Lc_fullsize(rect(1):rect(3), rect(2):rect(4)) = Lc;

    % get weighted center of cells in Lc_fullsize
    [fy, fx]= find(Lc_fullsize>0);
    center_cells_x = round( mean(fx) ); % before: round( (max(fx)+min(fx))/2 );
    center_cells_y = round( mean(fy) ); % before: round( (max(fy)+min(fy))/2 );

    % determine offset
    offset_x = round( size(Lc_fullsize,2)/2 - center_cells_x);
    offset_y = round( size(Lc_fullsize,1)/2 - center_cells_y);
    % write Lc into full size image
      % first, check if coordinates are out of bounds -> don't rescale image
      minrow=min(fy)+offset_y;
      maxrow=max(fy)+offset_y;
      mincol=min(fx)+offset_x;
      maxcol=max(fx)+offset_x;
      sizeLcfull=size(Lc_fullsize);
      if minrow<1 | maxrow>sizeLcfull(1) | mincol<1 | maxcol>sizeLcfull(2)
          fprintf('Image too large for centering (idx out of bounds). Will use non-centered image.')
          Lc_fullsize_centered=Lc_fullsize;
      else
          Lc_fullsize_centered( minrow:maxrow, mincol:maxcol ) = Lc_fullsize( min(fy):max(fy), min(fx):max(fx) );
      end

end
