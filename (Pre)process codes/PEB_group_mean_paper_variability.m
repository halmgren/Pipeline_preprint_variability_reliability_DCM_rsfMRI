function PEB_group_mean_paper_variability(SPM_dir,Work_dir)

all_ROI_defs={'Smith'};
all_procedure_names={'Basic','GSR','ROI_Size'};

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
        
        if strcmp(procedure,'Basic')
            mkdir([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
            
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
                                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_A_mean_' dataset '_subject_' num2str(subject) '_' ntwrk_name{network_number} '.mat'],'ME');
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
                        end
                        
                        clear PEB;
                    end
                end
                
                M = struct();
                M.alpha = 1;
                M.beta  = 16;
                M.hE    = 0;
                M.hC    = 1/16;
                M.Q     = 'all';
                M.X=ones(length(PEB_group),1);
                
                [PEB DCM]=spm_dcm_peb_of_peb(PEB_group',M,'A');
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group');
                
                clear PEB_group;
            end
        end
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')

            %%%%%%%%%%
            %Compare GSR basic
            %%%%%%%%%%
            mkdir([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
            for network_number=1:length(ntwrk_name)
                
                tmp=0;
                
                disp(network_number);
                for number_dataset=1:4
                    
                    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
                    for subject=1:number_subject
                        clear PEB;
                        tmp=tmp+1;
                        
                        %If no suprathresholod voxels!
                        try
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' ntwrk_name{network_number} '.mat']);
                        catch ME
                            tmp=tmp-1;
                            disp('PEB group')
                            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                             pause;
                            if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                                continue;
                            else
                                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_A_mean_comp_GSR_' dataset '_subject_' num2str(subject) '_' ntwrk_name{network_number} '.mat'],'ME');
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
                        end
                        
                        clear PEB;
                    end
                end
                
                M = struct();
                M.alpha = 1;
                M.beta  = 16;
                M.hE    = 0;
                M.hC    = 1/16;
                M.Q     = 'all';
                M.X=ones(length(PEB_group),1);
                
                [PEB DCM]=spm_dcm_peb_of_peb(PEB_group',M,'A');
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group')
                clear PEB_group;
                
            end
        end
    
        
        if strcmp(procedure,'ROI_Size')
            for size_ROI=1:4
                mkdir([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_group/']);
                for network_number=1:length(ntwrk_name)
                    tmp=0;
                    
                    disp(network_number);
                    for number_dataset=1:4
                        
                        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
                        for subject=1:number_subject
                            tmp=tmp+1;
                            try
                                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat']);
                            catch ME
                                tmp=tmp-1;
                                disp('PEB group')
                                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                                 pause;
                                if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                                    continue;
                                else
                                    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_group/ERROR_PEB_A_mean_comp_ROIsize_' dataset '_subject_' num2str(subject) '_' ntwrk_name{network_number} '.mat'],'ME');
                                    continue;
                                end
                            end
                            if length(PEB.Snames)<8
                                tmp=tmp-1;
                                continue
                            else
                                PEB_group{tmp}=PEB;
                            end
                            clear PEB;
                        end
                    end
                    
                    M = struct();
                    M.alpha = 1;
                    M.beta  = 16;
                    M.hE    = 0;
                    M.hC    = 1/16;
                    M.Q     = 'all';
                    M.X=ones(length(PEB_group),1);
                    
                    [PEB DCM]=spm_dcm_peb_of_peb(PEB_group',M,'A');
                    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group')
                    clear PEB_group;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%
        %Lateralization
        %%%%%%%%%%%%%%%%%
        if strcmp(procedure,'Basic')
            for network_number=1:length(ntwrk_name)
                load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB')
                
                %DMN
                if strcmp(ntwrk_name{network_number},'DMN')
                    T = 0;
                    C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
                    c = C'*spm_vec(PEB.Ep);
                    v = C'*PEB.Cp*C;
                    PP   = 1-spm_Ncdf(T,c,v);
                end
                
                mean_diff=c;
                var_of_sum=v;
                posterior_probability=PP;
                
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_A_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
            end
        end
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
            for network_number=1:length(ntwrk_name)
                load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'PEB')
                
                %DMN
                if strcmp(ntwrk_name{network_number},'DMN')
                    T = 0;
                    C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
                    c = C'*spm_vec(PEB.Ep);
                    v = C'*PEB.Cp*C;
                    PP   = 1-spm_Ncdf(T,c,v);
                end
                
                mean_diff=c;
                var_of_sum=v;
                posterior_probability=PP;
                
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_A_comp_GSR_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
            end
        end
        
        if strcmp(procedure,'ROI_Size')
            
            for network_number=1:length(ntwrk_name)
                for size_ROI=1:4
                    load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_group/PEB_A_mean_comp_ROIsize_group_' ntwrk_name{network_number} '.mat'],'PEB')
                   
                    %DMN
                    if strcmp(ntwrk_name{network_number},'DMN')
                        T = 0;
                        C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
                        c = C'*spm_vec(PEB.Ep);
                        v = C'*PEB.Cp*C;
                        PP   = 1-spm_Ncdf(T,c,v);
                    end
                    
                    mean_diff=c;
                    var_of_sum=v;
                    posterior_probability=PP;
                    
                    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_group/Lateralization_index_A_comp_ROIsize_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                    clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                end
            end
        end
    end
    
end
