function [timesPerframe_nrs, paramPerframe_nrs, schnitzesPerFrame_nrs, plottingfXaxis] = MW_getparamfromschnitzcells(schnitzcells, timeField, fieldOfInterest)
%function [timesPerframe_nrs, paramPerframe_nrs, schnitzesPerFrame_nrs, plottingfXaxis] = MW_getparamfromschnitzcells(schnitzcells, timeField, fieldOfInterest)
%
% Function that returns certain parameter value for each schnitz per frame
% based on MW_growthplots.m.
%
% Input:
% schnitzcells      standard schnitzcells struct
% timeField         name of the field to use for time parameter
% fieldOfInterest   idem for field of interest
%
% Output:
% - timesPerframe_nrs       gives applicable time for this frame
% - paramPerframe_nrs       gives values of parameter for this frame
% - schnitzesPerFrame       gives schnitz nrs for this frame (convenient 
%                           for color plotitng)
% - plottingfXaxis          array equal size of other output, but filled
%                           with applicable frame numbers, convenient for 
%                           plotting.
% Note on output:
% - note that the index for the output corresponds to frame numbers as
%                           given by 
%                           myframe_nrs = unique([schnitzcells.frame_nrs])
% 
% MW 2016.01

%% Obtain ranges from schnitzcells struct

myframe_nrs = unique([schnitzcells.frame_nrs]);
myTimes = unique([schnitzcells.(timeField)]);
allLengths = [schnitzcells.(fieldOfInterest)];

%%

% Surely this can be done prettier, but let's just loop to get the desired
% data now.
paramPerframe_nrs = {}; plottingfXaxis = {}; timesPerframe_nrs = {};
schnitzesPerFrame_nrs = {};
for f = myframe_nrs
   
    valueForThisFrame = [];
    timesForThisFrame = [];
    schnitzesForThisFrame = [];
    
    for schnitzIdx = 1:numel(schnitzcells)
        % Look whether this schnitz lives in current frame number,
        % and what timepoint belongs to the current frame.
        pointInSchnitz = find(schnitzcells(schnitzIdx).frame_nrs==f);
        
        % If so, add length at that timepoint to the collection of lengths
        if ~isempty(pointInSchnitz)
            timesForThisFrame(end+1) = schnitzcells(schnitzIdx).(timeField)(pointInSchnitz);
            valueForThisFrame(end+1) = schnitzcells(schnitzIdx).(fieldOfInterest)(pointInSchnitz);
            schnitzesForThisFrame(end+1) = schnitzIdx;
        end
    end
    
    timesPerframe_nrs{end+1} = timesForThisFrame;
    paramPerframe_nrs{end+1} = valueForThisFrame;
    schnitzesPerFrame_nrs{end+1} = schnitzesForThisFrame;
    plottingfXaxis{end+1} = ones(numel(valueForThisFrame),1)*f; % convenient as x-axis later
    
end

end