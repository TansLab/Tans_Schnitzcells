function slookup = MW_makeslookup(p)
% function slookup = MW_makeslookup(p)
% Code stolen from NW_makeMovie_branchGroups
%
% Returns lookup table slookup(frame, cellno) to obtain schnitznr.
%
% Note that the highest number of the label indices used corresponds to 
% maxCellNo = size(slookup,2);

%% load schnitzcells file
myFileName = [p.tracksDir p.movieName '_lin.mat'];
disp(['Loading ' myFileName '.']);
myData = load(myFileName);
s = myData.schnitzcells;

%% Create lookup table

% total number of individual cells (schnitzes)
nrOfSchnitzes = length(s);

% Loop over each individual cell
for i = 1:nrOfSchnitzes
  % Loop over the frames where it exists
  for j = 1:length(s(i).frame_nrs)
      
    if ~isempty(s(i).frame_nrs) 
        
      if s(i).frame_nrs
        
        % slookup(framenr,cellno)=schnitznr
        slookup(s(i).frame_nrs(j),s(i).cellno(j))=i;
        
      end
      
    end
    
  end
end

disp('Lookup table complete');

end