function Compute_BS_sign_stability_paper_variability(SPM_dir,Work_dir)

name_ROI_def='Smith';
[ROI_list]=Define_ROIs_paper_variability(name_ROI_def);
procedure='Basic';
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

for network_number=1:length(ntwrk_name)
    clear A_matrix;
    tmp=0;
    for number_dataset=1:4
        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
        for subject=1:number_subject
            
            clear PEB;
            tmp=tmp+1;
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
            catch
                continue;
            end
            %We don't include participants with less then 8
            %useful sessions.
            if length(PEB.Snames)<8
                tmp=tmp-1;
                continue
            else
                PEB_group{tmp}=full(vec2mat(PEB.Ep(1:16),4)');
                Lat_group{tmp}=mean_diff;
                ci=spm_invNcdf(1-0.05);
                EP=full(vec2mat(PEB.Ep(1:16),4)');
                CP=diag(PEB.Cp);
                CP=full(vec2mat(CP(1:16),4)');
                sgn_group(:,:,tmp)=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));
                
                PEB_group{tmp}(sgn_group(:,:,tmp)==-1)=NaN;
            
            end
            
            
            clear PEB;
            
        end
    end
    
    for subject=1:length(PEB_group)
        A_matrix(:,:,subject)=full(PEB_group{subject});
    end
    
    %Reorder matrices
    for subject=1:size(A_matrix,3)
        A_matrix_delat(:,:,subject)=A_matrix(:,:,subject);
    end
    
    for subject=1:length(Lat_group)
        if Lat_group{subject}>0
            A_matrix_delat([2 4],:,subject)=A_matrix_delat([4 2],:,subject);
            A_matrix_delat(:,[2 4],subject)=A_matrix_delat(:,[4 2],subject);
        end
    end
    
    for sz1=1:4
        for sz2=1:4
            number_nega(sz1,sz2)=sum(squeeze(sign(A_matrix_delat(sz1,sz2,:)))==-1);
            number_posi(sz1,sz2)=sum(squeeze(sign(A_matrix_delat(sz1,sz2,:)))==1);
            number_NaN(sz1,sz2)=sum(isnan(squeeze(A_matrix_delat(sz1,sz2,:))));
            
        end
    end
    
    save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BS_sign_stability_delat_' ntwrk_name{network_number} '.mat'],'number_nega','number_posi','number_NaN');
end

end