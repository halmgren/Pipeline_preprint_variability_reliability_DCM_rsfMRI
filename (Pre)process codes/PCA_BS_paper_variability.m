function PCA_BS_paper_variability(SPM_dir,Work_dir)

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
                    
                    try
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat']);
                    catch ME
                        tmp=tmp-1;
                        disp('PEB group')
                        disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                        %                             pause;
                        if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                            continue;
                        else
                            continue;
                        end
                    end
                    %We don't include participants with less then 8
                    %useful sessions.
                    if length(PEB.Snames)<8
                        tmp=tmp-1;
                        continue
                    else
                        PEB_group{tmp}=PEB;
                        Connectivity_group(:,tmp)=full(PEB.Ep);
                    end
                    
                    clear PEB;
                end
            end
            
%             clear PEB_group;
        end
        
    end
end

[coeff,score,latent,tsquared,explained,mu]=pca(Connectivity_group');
save([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/PCA_group.mat'],'coeff','score','latent','tsquared','explained','mu');

end