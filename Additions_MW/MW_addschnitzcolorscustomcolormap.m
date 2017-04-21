function p = MW_addschnitzcolorscustomcolormap(p)


highestSchnitzIndx = max(p.slookup(size(p.slookup,1),:));
% Set up custom colormap
% easy way:
%theColorMap = linspecer(maxCellNo);

% let's use one of the standard colormaps
standardColorMap = hsv(highestSchnitzIndx); % hsv jet

% but mix it up such that neighbors have different
% colors                    
shuffle=randperm(highestSchnitzIndx); 
%shuffle = [10  20  28  35  29  27  64  40  33  37  47  58  22  36  18  21  50  57  34  25  19  43   1  49  16  60  23   3  48   9  45  38  44  46   4  26   7  15  54  59  55  24   6  14  42  13   8  61  63  41   5   2  30  53  31  52  32  51  56  17  11  62  12  39];

% perform shuffling
standardColorMapShuffled = NaN(size(standardColorMap));
standardColorMapShuffled(:,1) = standardColorMap(shuffle,1);
standardColorMapShuffled(:,2) = standardColorMap(shuffle,2);
standardColorMapShuffled(:,3) = standardColorMap(shuffle,3);

% how many copies of the colormap do we need?
%copiesNeeded = ceil(maxCellNo/COLORMAPSIZE);

% create the color map
%theColorMap = repmat(standardColorMapShuffled,copiesNeeded,1);

% set the colormap
p.customColors = [0 0 0; standardColorMapShuffled; 1 1 1];

end