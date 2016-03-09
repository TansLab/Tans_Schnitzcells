function Iout = PN_reseed(Iin,phim,ppx,ppy)

Iout = Iin;



cIout = imcomplement(Iout);
se = strel('disk',1);
IoutErode = imerode(cIout, se);                                       %Morphological reconstruction ...
cIout = imreconstruct(IoutErode,cIout);                             %... in the mask of the original image
cIout = imdilate(cIout, se);                                           %Allows a better edge determination at contacts.

edgeIm = edge(cIout,'log',0,2);

newCell = imfill(edgeIm,[ppx ppy]) & ~edgeIm;
newCell = imdilate(newCell,strel('diamond',2));

newCell = newCell & ~(logical(Iout));

Iout(newCell) = max2(Iin)+1;



end