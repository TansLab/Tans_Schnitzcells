function [Min XMin YMin] = MinCoordinates2_N(matr,N)
% same functionality as MinCoordinates2 but returns coordinates vectors with the N lowest
% distances
% if e.g. N==1 and 4 points have the same value -> a vector of size 4 is
% returned
%
% matr has to be >=0

idxneg=find(matr(:)<0);
if ~isempty(idxneg)
    error('no negative entries in matrix allowed.')
end

matrneg=-matr;
allvalues=unique(matrneg(:));
choosevalues=allvalues(max(1,end-N+1):end);

[idxrow idxcol]=find(ismember(matrneg,choosevalues));
idxlinear=find(ismember(matrneg,choosevalues));


Min=matr(idxlinear);
[Min,idxorder]=sort(Min);

XMin=idxcol(idxorder);
YMin=idxrow(idxorder);