
% AutoCorr & Xcorr with errorbars (Here, noise is calculated seperately for each branchgroup in the array of branchgroups!)
fitTime = fitTime + [2 -2];
branchData2 = DJK_getBranches(p,s_rm,'dataFields',{'R_time'    'R6_mean_cycCor' ,'G6_mean_cycCor' 'muP15_fitNew_cycCor' }, 'fitTime', fitTime); 
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_NoisePerBranchGroup'];
trimmed_branches2 = DJK_trim_branch_data(branchData2);
branch_groups2 = DJK_divide_branch_data(trimmed_branches2);
branch_groups2_noise = NW_addToBranchGroups_noise(p, branch_groups2,'dataFields',{'R_time' 'R6_mean_cycCor'  'G6_mean_cycCor' 'muP15_fitNew_cycCor'});
branch_groups=branch_groups2_noise;

DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_G6_mean_cycCor', 'noise_muP15_fitNew_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_R6_mean_cycCor', 'noise_muP15_fitNew_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_G6_mean_cycCor',  'noise_R6_mean_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_G6_mean_cycCor',  'noise_G6_mean_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_R6_mean_cycCor',  'noise_R6_mean_cycCor','selectionName',name_rm_branch,'timeField','R_time');



branchData2 = DJK_getBranches(p,s_rm,'dataFields',{'dR5_time'  'R_time'   'muP15_fitNew_atdR5_cycCor' 'dR5_cycCor'  'dG5_cycCor'  }, 'fitTime', fitTime); 
name_rm_branch = [name_rm '_' num2str(fitTime(1)) '_' num2str(fitTime(2)) '_NoisePerBranchGroup'];
trimmed_branches2 = DJK_trim_branch_data(branchData2);
branch_groups2 = DJK_divide_branch_data(trimmed_branches2);
branch_groups2_noise = NW_addToBranchGroups_noise(p, branch_groups2,'dataFields',{'dR5_time'  'R_time'   'muP15_fitNew_atdR5_cycCor' 'dR5_cycCor'  'dG5_cycCor' });
branch_groups=branch_groups2_noise;

DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dG5_cycCor', 'noise_muP15_fitNew_atdR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dR5_cycCor', 'noise_muP15_fitNew_atdR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dG5_cycCor',  'noise_dR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dG5_cycCor',  'noise_dG5_cycCor','selectionName',name_rm_branch,'timeField','R_time');
DJK_plot_crosscorrelation_standard_error_store(p, branch_groups, 'noise_dR5_cycCor',  'noise_dR5_cycCor','selectionName',name_rm_branch,'timeField','R_time');

