function [timesPerframe_nrs, lengthsPerframe_nrs, plottingfXaxis] = MW_growthfromschnitzcells(schnitzcells, timeField, lengthField)
% function [timesPerframe_nrs lengthsPerframe_nrs plottingfXaxis] = MW_growthfromschnitzcells(schnitzcells, timeField, lengthField)
%
% Function that returns lengths for each schnitz per frame
% based on MW_growthplots.m.
%
% Input:
% schnitzcells      standard schnitzcells struct
% timeField         name of the field to use for time parameter
% lengthField       idem for length
%
% Output:
% timesPerframe_nrs
% lengthsPerframe_nrs
% plottingfXaxis
% 
% MW 2016.01

%% Obtain ranges from schnitzcells struct

myframe_nrs = unique([schnitzcells.frame_nrs]);
myTimes = unique([schnitzcells.(timeField)]);
allLengths = [schnitzcells.(lengthField)];

%%

% Surely this can be done prettier, but let's just loop to get the desired
% data now.
lengthsPerframe_nrs = {}; plottingfXaxis = {}; timesPerframe_nrs = {};
for f = myframe_nrs
   
    lengthForThisFrame = [];
    timesForThisFrame = [];
    
    for schnitzIdx = 1:numel(schnitzcells)
        % Look whether this schnitz lives in current frame number,
        % and what timepoint belongs to the current frame.
        pointInSchnitz = find(schnitzcells(schnitzIdx).frame_nrs==f);
        
        % If so, add length at that timepoint to the collection of lengths
        if ~isempty(pointInSchnitz)
            timesForThisFrame(end+1) = schnitzcells(schnitzIdx).(timeField)(pointInSchnitz);            
            lengthForThisFrame(end+1) = schnitzcells(schnitzIdx).(lengthField)(pointInSchnitz);
        end
    end
    
    timesPerframe_nrs{end+1} = timesForThisFrame;
    lengthsPerframe_nrs{end+1} = lengthForThisFrame;
    plottingfXaxis{end+1} = ones(numel(lengthForThisFrame),1)*f; % convenient as x-axis later
    
end

end