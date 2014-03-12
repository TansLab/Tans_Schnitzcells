function imnew = NW_renumberimage(imold);
% faster than original Version because only renumbers cells whose index is
% too high (NW 2012-05-10)

imnew = imold;
u = unique(imold(:));

%for i = 2:length(u), % loop costs 0.49 sec (for 400 cells)
%    imnew(imold==u(i))=i-1;
%end;


% cell numbers which are too high (=higher than #cells)
exceedNumbers=u(u>(length(u)-1))';

if ~isempty(exceedNumbers) % have to renumber
    % perfect numbering would be [0,1,2,3,4,5.... # cells]
    perfectnumbers=[0:1:length(u)];
    % find not used numbers which will be associated to cells with too high
    % cell number
    usedNumbers=intersect(u,perfectnumbers);
    idx=ismember(perfectnumbers,usedNumbers);
    freeNumbers=perfectnumbers(idx==0);
    
    runner=1;
    for i=exceedNumbers
        imnew(imold==i)=freeNumbers(runner);
        runner=runner+1;
    end
end