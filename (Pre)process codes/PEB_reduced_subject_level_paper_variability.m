function PEB_reduced_subject_level_paper_variability(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute hemispheric (or front-back) dominance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp=0;
name_ROI_def='Smith';
procedure='BMR';

[ROI_list]=Define_ROIs_paper_variability(name_ROI_def);

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

%%%%%%%%%%%%%%%%
%Lateralization
%%%%%%%%%%%%%%%%
for network_number=1:length(ntwrk_name)
    disp(ntwrk_name{network_number});
    
    mkdir([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
     
    load([Work_dir '/Results_paper_variability/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
    
    for subject=1:length(DCM)
                    
            %DMN
            if(strcmp(ntwrk_name{network_number},'DMN'))
                T = 0;
                C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(DCM{subject}.Ep)-16)]'/3;
                c(subject) = C'*spm_vec(DCM{subject}.Ep);
                v(subject) = C'*DCM{subject}.Cp*C;
                PP(subject)   = 1-spm_Ncdf(T,c(subject),v(subject));
            end
            
    end
    mean_diff=c;
    var_of_sum=v;
    posterior_probability=PP;
    
    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
    clear mean_diff var_of_sum c v C PP T posterior_probability;
    clear DCM PEB;
end

%%%%%%%%%%%%%%%%%
%Correlation
%%%%%%%%%%%%%%%%%
for network_number=1:length(ntwrk_name)
    disp(ntwrk_name{network_number});
    
    load([Work_dir '/Results_paper_variability/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
    
    sz=sqrt(length(DCM{1}.Ep));
    A_matrix=zeros(4,4,length(DCM));
    for subject=1:length(DCM)
        A_matrix(:,:,subject)=vec2mat(full(DCM{subject}.Ep),sz)';
    end
    
    B=nchoosek(1:length(DCM),2);
    for combi=1:length(B)
        tmp3=corrcoef(reshape(A_matrix(:,:,B(combi,1)),sz^2,1),reshape(A_matrix(:,:,B(combi,2)),sz^2,1));
        BS_correlation_BMR(combi)=tmp3(1,2);
    end
    
    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BS_correlation_' ntwrk_name{network_number} '.mat'],'BS_correlation_BMR');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correlation reordered matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for network_number=1:length(ntwrk_name)
    disp(ntwrk_name{network_number});
    
    load([Work_dir '/Results_paper_variability/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
    load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_' ntwrk_name{network_number} '.mat']);
    
    sz=sqrt(length(DCM{1}.Ep));
    A_matrix=zeros(4,4,length(DCM));
    for subject=1:length(DCM)
        A_matrix(:,:,subject)=vec2mat(full(DCM{subject}.Ep),sz)';
    end
    
    A_matrix_delat=A_matrix;
    
    for subject=1:length(DCM)
        if mean_diff(subject)>0
            A_matrix_delat([2 4],:,subject)=A_matrix_delat([4 2],:,subject);
            A_matrix_delat(:,[2 4],subject)=A_matrix_delat(:,[4 2],subject);
        end
    end
    
    B=nchoosek(1:length(DCM),2);
    for combi=1:length(B)
        tmp3=corrcoef(reshape(A_matrix_delat(:,:,B(combi,1)),[1 16]),reshape(A_matrix_delat(:,:,B(combi,2)),[1 16]));
        BS_correlation_delat_BMR(combi)=tmp3(1,2);
    end
    
    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BS_correlation_delat_' ntwrk_name{network_number} '.mat'],'BS_correlation_delat_BMR');
end

end


