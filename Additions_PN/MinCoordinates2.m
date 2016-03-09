function [Min XMin YMin] = MinCoordinates2(matr)

[values,Ypositions] = min(matr);
[Min,XMin] = min(values);
YMin = Ypositions(XMin);
