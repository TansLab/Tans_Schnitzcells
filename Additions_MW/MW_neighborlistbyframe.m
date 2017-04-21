
% MW, 2-3-2017

% Currently based on cellno's and segfile only

load('H:\EXPERIMENTAL_DATA_2017\2017-02-02_OAA_pulsing_3\pos2crop\segmentation\pos2cropseg040.mat')

% Cutoff distance in pixels
EXPANDCELLSBY =11;

% Select cellnumbers in this frame (note you need to convert them to schnitz
% nrs later if desired)
allCellNumbersThisFrame=unique(Lc(:));
allCellNumbersThisFrame=allCellNumbersThisFrame(allCellNumbersThisFrame>0);
% Create neighbor list
neighborlist={};
for cellIdx=allCellNumbersThisFrame'

    % Select current cell
	selectedCellMatrix=(Lc==cellIdx);
    % Dilate that cell
	se = strel('disk',EXPANDCELLSBY);
	selectedCellMatrixExpanded = imdilate(selectedCellMatrix,se);

    % Find out which neighbors that cell has (note: includes itself)
	neighborMatrix= selectedCellMatrixExpanded.*Lc;
    % Add these neighbors to a neighbor list
	neighborlist{cellIdx}=unique(neighborMatrix(:));

end

% neighborlist{i} now gives all neighbors for cell i