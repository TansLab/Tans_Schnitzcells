function slookup = MW_makeslookup(p)
% function slookup = MW_makeslookup(p)
% Code stolen from NW_makeMovie_branchGroups
%
% Returns lookup table slookup(frame, cellno) to obtain 

%% load schnitzcells file
myFileName = [p.tracksDir p.movieName '_lin.mat'];
disp(['Loading ' myFileName '.']);
myData = load(myFileName);
s = myData.schnitzcells;

%% Create lookup table

for i = 1:length(s)
  for j = 1:length(s(i).frame_nrs)
    if ~isempty(s(i).frame_nrs)
      if s(i).frame_nrs
        slookup(s(i).frame_nrs(j),s(i).cellno(j))=i;
      end
    end
  end
end

disp('Lookup table complete');

end