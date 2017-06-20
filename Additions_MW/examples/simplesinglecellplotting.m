

SCHNITZNR=45;

cellDataTime   = [schnitzcells(SCHNITZNR).time]
cellDataGrowth5 = [schnitzcells(SCHNITZNR).muP5_fitNew_all]
cellDataGrowth9 = [schnitzcells(SCHNITZNR).muP9_fitNew_all]
cellDataLength = [schnitzcells(SCHNITZNR).length_fitNew]
cellDataLengthSkeleton = [schnitzcells(SCHNITZNR).length_skeleton]

figure(1); clf; hold on;
plot(cellDataTime,cellDataLength,'o-','LineWidth',2);
plot(cellDataTime,cellDataLengthSkeleton,'o-','LineWidth',2);

legend('Polynomial fit','Skeleton length');
xlabel('Time (min)');
ylabel('Length');
MW_makeplotlookbetter(20);

figure(2); clf; hold on;
plot(cellDataTime,cellDataGrowth5,'-o','LineWidth',2);
plot(cellDataTime,cellDataGrowth9,'-o','LineWidth',2);

ylim([0,2]);

legend({'5 frames exp. fit','9 frames fit'});
xlabel('Time (min)');
ylabel('Growth rate (dbl/hr)');
MW_makeplotlookbetter(20);