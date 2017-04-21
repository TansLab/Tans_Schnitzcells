

% NOTE: I don't think this works..

FIELDOFINTEREST='length_fitNew';

myFrameIdxs = ourSettings.currentFrameRange;

%% gather data

yValuesArray = [];
for fridx = myFrameIdxs
    
    indicesThisFrame = ([schnitzcells.frame_nrs]==fridx);
    yValueSummed = sum([schnitzcells(indicesThisFrame).(FIELDOFINTEREST)]); 
   
    yValuesArray(end+1) = yValueSummed;
        
end

disp('done')

%% plot

figure(1); clf; hold on;
plot(myFrameIdxs , yValuesArray,'o')
set(gca,'yscale','log')
disp('done')





