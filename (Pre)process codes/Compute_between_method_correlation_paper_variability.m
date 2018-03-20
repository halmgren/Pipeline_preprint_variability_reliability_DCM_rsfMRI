function Compute_between_method_correlation_paper_variability(SPM_dir,Work_dir)

name_ROI_def='Smith';
[ROI_list]=Define_ROIs_paper_variability(name_ROI_def);

tmp=0;

for VOI_number=1:size(ROI_list,1)
    ntwrk=ROI_list{VOI_number,1}(1:3);
    
    if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
        ntwrk_size(tmp)=ntwrk_size(tmp)+1;
        continue
        
    else
        tmp=tmp+1;
        ntwrk_size(tmp)=1;
        ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
    end
end

%%%%%%%%%%%%%%%%%
%group level PEB
%%%%%%%%%%%%%%%%%
for network_number=1:length(ntwrk_name)
    comb_wo_rplc=nchoosek(1:4,2);
    multiproc_PEB{1}=load([Work_dir '/Results_paper_variability/DCM/ROI_Size/' name_ROI_def '/Full_model/4/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB');
    multiproc_PEB{2}=load([Work_dir '/Results_paper_variability/DCM/ROI_Size/' name_ROI_def '/Full_model/8/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB');
    multiproc_PEB{3}=load([Work_dir '/Results_paper_variability/DCM/ROI_Size/' name_ROI_def '/Full_model/12/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB');
    multiproc_PEB{4}=load([Work_dir '/Results_paper_variability/DCM/ROI_Size/' name_ROI_def '/Full_model/16/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB');
    
    for combi=1:size(comb_wo_rplc,1)
        correl=corrcoef(multiproc_PEB{comb_wo_rplc(combi,1)}.PEB.Ep(1:ntwrk_size(network_number)^2),multiproc_PEB{comb_wo_rplc(combi,2)}.PEB.Ep(1:ntwrk_size(network_number)^2));
        BM_correlation(combi)=correl(1,2);
    end
    
    save([Work_dir '/Results_paper_variability/DCM/ROI_Size/' name_ROI_def '/Full_model/BM_correlation.mat'],'BM_correlation');
    
    clear BM_correlation multiproc_PEB;
    
    %GSR
    comb_wo_rplc=nchoosek(1:2,2);
    multiproc_PEB{1}=load([Work_dir '/Results_paper_variability/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'PEB');
    multiproc_PEB{2}=load([Work_dir '/Results_paper_variability/DCM/GSR/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'PEB');
    
    for combi=1:size(comb_wo_rplc,1)
        correl=corrcoef(multiproc_PEB{comb_wo_rplc(combi,1)}.PEB.Ep(1:ntwrk_size(network_number)^2),multiproc_PEB{comb_wo_rplc(combi,2)}.PEB.Ep(1:ntwrk_size(network_number)^2));
        BM_correlation(combi)=correl(1,2);
    end
     
    save([Work_dir '/Results_paper_variability/DCM/GSR/' name_ROI_def '/Full_model/PEB_group/BM_correlation.mat'],'BM_correlation');
    
    clear BM_correlation multiproc_PEB;
    
end
