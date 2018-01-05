if ~isfield(p,'sameLength')
    p.sameLength=1; % assume branches are of equal length
end

dx = CorrData(2,1)-CorrData(1,1);
surfaceUnderCorrelation = sum(CorrData(:,2))*dx % note to self: CorrData is already only right half if autocorr
tau = surfaceUnderCorrelation;

mytimepoints = CorrData(:,1);
fittedExponential = exp(mytimepoints./-tau);

myxlim=[min(CorrData(:,1)), max(CorrData(:,1))];

h20=figure('Visible',FIGUREVISIBLE), clf, hold on
plot(myxlim, [0 0],'k-');% axis at zero
if p.sameLength
    errorbar(CorrData(:,1),CorrData(:,2),CorrData(:,3),'o','LineWidth',3,'Color',[.5 .5 .5])
end
plot(CorrData(:,1),CorrData(:,2),'o','LineWidth',3,'Color','k')
plot(mytimepoints,fittedExponential,'-','LineWidth',3,'Color',[.2 .2 1])

xlim(myxlim);
title([p.movieDate ', ' p.movieName 10 'Area under curve, \tau: ' num2str(surfaceUnderCorrelation)])
xlabel('\tau [hr]')
ylabel([associatedFieldNames{2} 10 'autocorrelation [normalized]'],'Interpreter','None');
MW_makeplotlookbetter(20)

saveas(h20,[myOutputFolder 'TIF_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_expFit.tif']);
saveas(h20,[myOutputFolder 'EPS_crosscorrs_' associatedFieldNames{1,2} '_' associatedFieldNames{1,3} '_expFit.eps'],'epsc');
