
function X = readphaseset(pname, totalnum)
% function X = readphaseset(pname, totalnum)
%
% reads set of phase images
%  pname:   'coliglo/2002-04-01/LacSymBranch-03-p-2-046.tif'
%  totalnum: total number of phase slices taken

f = findstr(pname,'-p-');
xtra = 0;

if totalnum==1
    for i = 1:3,
        name = pname;
        im = imread(name);
        X(:,:,i) = im;
    end
else
    for i = 1:totalnum,
        name = pname;
        name(f+3) = num2str(i+xtra);
        im = imread(name);
        X(:,:,i) = im;
    end
end;
