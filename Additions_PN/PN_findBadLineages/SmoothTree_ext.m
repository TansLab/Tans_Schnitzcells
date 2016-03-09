function schnitzcells = SmoothTree(schnitzcells,fieldName,span)
% !!!uses the function smooth of matlab
% NOT MODIFIED YET!!! DONT USE !!!
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
        parentsVec_ext = AscendField_ext(schnitzcells,fieldName,i);
        vectParents = parentVec_ext( max(1,end-halfspan+1):end );
        
        %forward extension of the vector
        kidsVec_ext = DescendAverage_ext(schnitzcells,fieldName,i);
        vectKids = kidsVec_ext( 1:min(length(kidsVec_ext),halfspan) );
  
        %create the extended then smoothed vector
        vectTot = [vectParents(:) ; svect(:) ; vectKids(:)]';
        smoothTemp = smooth(vectTot,span);
        iniPoint = length(vectParents);
        s.(newField) = (smoothTemp( iniPoint + (1:length(svect)) ))';
        schnitzcells(i) = s;
        
    end
end



end
