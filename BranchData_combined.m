% COMBINE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
branchDataCombined = struct;
for i = 1:length(branchData1)
  branchDataCombined(i).schnitzNrs = branchData1(i).schnitzNrs;
  branchDataCombined(i).Y_time = branchData1(i).Y_time;
  branchDataCombined(i).Y6_mean = branchData1(i).Y6_mean;
  branchDataCombined(i).weight = branchData1(i).weight;
  branchDataCombined(i).muP17_fitNew= branchData1(i).muP17_fitNew;
  branchDataCombined(i).noise_Y6_mean = branchData1(i).noise_Y6_mean;
  branchDataCombined(i).noise_muP17_fitNew= branchData1(i).noise_muP17_fitNew;
end
for i = 1:length(branchData2)
  branchDataCombined(length(branchData1)+i).schnitzNrs = branchData2(i).schnitzNrs;
  branchDataCombined(length(branchData1)+i).Y_time = branchData2(i).Y_time;
  branchDataCombined(length(branchData1)+i).Y6_mean = branchData2(i).Y6_mean;
  branchDataCombined(length(branchData1)+i).weight = branchData2(i).weight;
  branchDataCombined(length(branchData1)+i).muP17_fitNew= branchData2(i).muP17_fitNew;
  branchDataCombined(length(branchData1)+i).noise_Y6_mean = branchData2(i).noise_Y6_mean;
  branchDataCombined(length(branchData1)+i).noise_muP17_fitNew= branchData2(i).noise_muP17_fitNew;
end
