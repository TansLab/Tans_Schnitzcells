function schnitzcells = SmoothTree(schnitzcells,fieldName,span)
% !!!uses the function smooth of matlab

if mod(span,2) == 0
    error('The span should be an odd number.')
end
halfspan = int8((span-1)/2);
ls = length(schnitzcells);
newField = [fieldName '_smoothed'];

%initialization
for i = 1:ls
    schnitzcells(i).(newField) = [];
end

for i = 1:ls
    s = schnitzcells(i);
    svect = s.(fieldName);
   
    if ~isempty(svect)
        
        %backward extension of the vector
        parentsAverage = AscendField(schnitzcells,fieldName,i);
        vectParents = parentsAverage( max(1,end-halfspan+1):end );
        
        %forward extension of the vector
        kidsAverage = DescendAverage(schnitzcells,fieldName,i);
        vectKids = kidsAverage( 1:min(length(kidsAverage),halfspan) );
  
        %create the extended then smoothed vector
        vectTot = [vectParents(:) ; svect(:) ; vectKids(:)]';
        smoothTemp = smooth(vectTot,span);
        iniPoint = length(vectParents);
        s.(newField) = (smoothTemp( iniPoint + (1:length(svect)) ))';
        schnitzcells(i) = s;
        
    end
end



end
