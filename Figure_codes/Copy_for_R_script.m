function Copy_for_R_script(SPM_dir,Work_dir)

load([Work_dir '/DatasetKirby/sub-01_results/DCM/Basic/Smith/Full_model/PEB_A_mean_DMN.mat'],'PEB');
S1=full(reshape(PEB.Ep,4,4));
load([Work_dir '/DatasetKuehn/sub-01_results/DCM/Basic/Smith/Full_model/PEB_A_mean_DMN.mat'],'PEB');
S2=full(reshape(PEB.Ep,4,4));
load([Work_dir '/DatasetGordon/sub-06_results/DCM/Basic/Smith/Full_model/PEB_A_mean_DMN.mat'],'PEB');
S3=full(reshape(PEB.Ep,4,4));
load([Work_dir '/DatasetGordon/sub-07_results/DCM/Basic/Smith/Full_model/PEB_A_mean_DMN.mat'],'PEB');
S4=full(reshape(PEB.Ep,4,4));

save([Work_dir '/Figures_paper_variability/Figure2_matrices.mat'],'S1','S2','S3','S4');

end