function PCA_WS_paper_variability(SPM_dir,Work_dir)

all_procedure_names={'Basic'};
all_ROI_defs={'Smith'};

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROI_list]=Define_ROIs_paper_variability(name_ROI_def);
    
    tmp1=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            ntwrk_size(tmp1)=ntwrk_size(tmp1)+1;
            continue
            
        else
            tmp1=tmp1+1;
            ntwrk_size(tmp1)=1;
            ntwrk_name{tmp1}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1:length(all_procedure_names)
        procedure=all_procedure_names{number_procedure};
        
        
        for network_number=1:length(ntwrk_name)
            
            tmp=0;
            
            disp(network_number);
            for number_dataset=1:4
                
                [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
                for subject=1:number_subject
                    clear PEB;
                    tmp=tmp+1;
                    
                    if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
                        tmp=tmp-1;
                    continue
                end
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                    
                    flag=zeros(1,length(GCM));
                    for diagn=1:length(GCM)
                        if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                            flag(diagn)=1;
                        end
                    end
                    
                    GCM(find(flag==1))=[];
                    
                    GCM_group{tmp}=GCM;
                    for session=1:length(GCM)
                        Connectivity_group{tmp}(session,:)=GCM{session}.Ep.A(:);
                    end

                    
                    clear PEB;
                end
            end
            
        end
        
    end
end

for subject=1:length(Connectivity_group)
    [coeff{subject},score{subject},latent{subject},tsquared{subject},explained{subject},mu{subject}]=pca(Connectivity_group{subject});
end

save([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/PCA_subject_spec.mat'],'coeff','score','latent','tsquared','explained','mu');
end