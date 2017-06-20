
%% Some general stats
nrBarrenIssues = size(problemSchnitzesBarren,2);
nrMovingIssues =  size(problemSchnitzesMoving ,2);
nrSlowGrowthIssues = size(problemSchnitzesSlowGrowth  ,2);
nrFastGrowthIssues = size(problemSchnitzesFastGrowth  ,2);
nrDivisionIssues = size(problemSchnitzesDivLenChange  ,2);

disp('Some summary stats:');
disp(['nrBarrenIssues = ' num2str(nrBarrenIssues)]);
disp(['nrMovingIssues = ' num2str(nrMovingIssues)]);
disp(['nrSlowGrowthIssues = ' num2str(nrSlowGrowthIssues)]);
disp(['nrFastGrowthIssues = ' num2str(nrFastGrowthIssues)]);
disp(['nrDivisionIssues = ' num2str(nrDivisionIssues)]);

%% identify non-barren related suspicious cells
theBarrenOnes = ismember(p.problemCells(:,1)',problemSchnitzesBarren);
nonBarrenIssues = ~theBarrenOnes;

problemCells=p.problemCells;
problemCellSelection = problemCells(nonBarrenIssues,:);

%% update problemcells such that only those are kept

problemCellSelection = problemCells(nonBarrenIssues,:);
p.problemCellsOld=p.problemCells;
p.problemCells=problemCellSelection;

disp('Removed barren cells from problem set.');