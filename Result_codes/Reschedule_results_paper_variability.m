function Reschedule_results_paper_variability(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subject-spec PEB results and asymmetry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
                        catch ME
                            tmp=tmp-1;
                            disp('PEB group')
                            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                           
                        end
                        %We don't include participants with less then 8
                        %useful sessions.
                        if length(PEB.Snames)<8
                            tmp=tmp-1;
                            continue
                        else
                            PEB_group{tmp}=PEB;
                            mean_diff_group(tmp)=mean_diff;
                            var_of_sum_group(tmp)=var_of_sum;
                            posterior_probability_group(tmp)=posterior_probability;
                        end
                        
                        clear PEB;
                    end
                end
                
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Group_array_results_' ntwrk_name{network_number} '.mat'],'PEB_group','mean_diff_group','var_of_sum_group','posterior_probability_group');
                
                clear PEB_group;
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%sign-stability
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count positives and negatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure='Basic';
name_ROI_def='Smith';
ntwrk_name='DMN';

clear A_matrix;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name '.mat']);
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name '.mat']);
        catch ME
            tmp=tmp-1;
            disp('PEB group')
            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
            %
            continue;
            
        end
        %We don't include participants with less then 8
        %useful sessions.
        if length(PEB.Snames)<8
            tmp=tmp-1;
            continue
        else 
            
            PEB_group{tmp}=PEB;
            Lat_group{tmp}=mean_diff;
        end
        A_matrix_tmp=full(vec2mat(PEB.Ep(1:16),4)');
        
        %CI
        ci=spm_invNcdf(1-0.05);
        EP=full(vec2mat(PEB.Ep(1:16),4)');
        CP=diag(PEB.Cp);
        CP=full(vec2mat(CP(1:16),4)');
        sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));
        
        A_matrix_tmp(sgn==-1)=NaN;
        A_matrix(:,:,tmp)=A_matrix_tmp;
        clear ci EP CP sgn;
    end
end

A_matrix_delat=A_matrix;


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

%Make second column and row dominant and fourth non-dominant
number_nega([2 4],:)=number_nega([4 2],:);
number_nega(:,[2 4])=number_nega(:,[4 2]);
number_posi([2 4],:)=number_posi([4 2],:);
number_posi(:,[2 4])=number_posi(:,[4 2]);
number_NaN([2 4],:)=number_NaN([4 2],:);
number_NaN(:,[2 4])=number_NaN(:,[4 2]);

%Ignore self-connections
number_posi(logical(eye(size(number_posi))))=0;
number_nega(logical(eye(size(number_nega))))=0;
number_NaN(logical(eye(size(number_NaN))))=0;

%proportion iso counts
prop_posi=number_posi/tmp;
prop_nega=number_nega/tmp;
prop_NaN=number_NaN/tmp;

cmap=[1 1 1; 0 1/2 0; 0 1 0; 1/2 0 0; 1 0 0];

sign_stab=zeros(4,4);
for sz1=1:4
    for sz2=1:4
        if prop_posi(sz1,sz2)>0.8
            sign_stab(sz1,sz2)=1;
        elseif prop_posi(sz1,sz2)>0.6
            sign_stab(sz1,sz2)=2;
        elseif prop_nega(sz1,sz2)>0.8
            sign_stab(sz1,sz2)=3;
        elseif prop_nega(sz1,sz2)>0.6
            sign_stab(sz1,sz2)=4;
        end
    end
end

save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Sign_consistency.mat'],'number_nega','number_posi','number_NaN','prop_posi','prop_nega','prop_NaN','sign_stab','A_matrix_delat');

%%%%%%%%%%%%%%%%%%%%%%%%
%Reliability
%%%%%%%%%%%%%%%%%%%%%%%%

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/WS_correlation_DMN.mat']);
        cr_group(1,tmp)=mean(cr_6);
        
    end
end
cr_group(:,end+1)=mean(cr_group,2);

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/WS_correlation_delat_DMN.mat']);
        cr_group_delat(1,tmp)=mean(cr_6);
    end
end

cr_group_delat(:,end+1)=mean(cr_group_delat,2);

cr_group_both=[cr_group;cr_group_delat];

save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WS_correlation.mat'],'cr_group','cr_group_delat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reliability Asymmetry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mean_diff_group;

tmp=0;
tmp2=-1;

for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        tmp2=tmp2+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These specific subjects were excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            tmp2=tmp2-1;
            continue
        end
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/Lateral_index_individ_DMN.mat']); %load lateralization
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']); %load QC's
        
        
        flag=zeros(1,length(mean_diff));
        for diagn=1:length(mean_diff)
            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                flag(diagn)=1;
            end
            
        end
        
        mean_diff(flag==1)=[];
        posterior_probability(flag==1)=[];
        flag2=zeros(1,length(mean_diff));
        N_sessions(tmp)=length(posterior_probability);
        
        for diagn=1:length(mean_diff)
            if posterior_probability(diagn)<0.95&&posterior_probability(diagn)>0.05
                flag2(diagn)=1;
            end
        end
        
        mean_diff_group{tmp}=mean_diff;
        stab(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff);
        Non_sign{tmp,:}=flag2;
    end
end

[~,order]=sort(stab);
for i=1:length(stab)
    Non_sign_order{i}=Non_sign{order(i)};
    mean_diff_group_order{i}=mean_diff_group{order(i)};
    
    stab2(i)=stab(order(i));
end

save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Asymmetry_reliability.mat'],'stab2','stab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WS cross-correlation BMR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/BMR/Smith/Full_model/WS_correlation_DMN.mat']);
        cr_group_BMR(1,tmp)=mean(cr_6);
        
    end
end
cr_group_BMR(:,end+1)=mean(cr_group_BMR,2);

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/BMR/Smith/Full_model/WS_correlation_delat_DMN.mat']);
        cr_group_delat_BMR(1,tmp)=mean(cr_6);
    end
end

cr_group_delat_BMR(:,end+1)=mean(cr_group_delat_BMR,2);

cr_group_both=[cr_group_BMR;cr_group_delat_BMR];

save([Work_dir '/Results_paper_variability/DCM/BMR/Smith/Full_model/PEB_group/WS_correlation_BMR.mat'],'cr_group_BMR','cr_group_delat_BMR');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asymmetry: subject-level
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mean_diff_group;
clear var_of_sum_group;
clear posterior_probability_group;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        for size_ROI=1:4
            
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/ROI_Size/Smith/Full_model/' num2str(size_ROI*4) '/Lateralization_index_A_comp_ROIsize_DMN.mat']);
            mean_diff_group(size_ROI,tmp)=mean_diff;
            var_of_sum_group(size_ROI,tmp)=var_of_sum;
            posterior_probability_group(size_ROI,tmp)=posterior_probability;
        end
    end
end

save([Work_dir '/Results_paper_variability/DCM/ROI_Size/Smith/Full_model/Subject_asymmetry.mat'],'mean_diff_group','var_of_sum_group','posterior_probability_group');

%%%%%%%%%%%%%%%%%%%%%%%
%Reliability: ROI size
%%%%%%%%%%%%%%%%%%%%%%%
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        for size_ROI=1:4
            
            
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/ROI_Size/Smith/Full_model/' num2str(size_ROI*4) '/WS_correlation_comp_ROIsize_DMN.mat']);
            cr_group_ROIsize(size_ROI,tmp)=mean(cr_6);
        end
    end
end

cr_group_ROIsize(:,end+1)=mean(cr_group_ROIsize,2);

save([Work_dir '/Results_paper_variability/DCM/ROI_Size/Smith/Full_model/WS_correlation_ROIsize.mat'],'cr_group_ROIsize');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%GSR: Asymmetry: subject-level
%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/Lateralization_index_A_comp_GSR_DMN.mat']);
        mean_diff_group(1,tmp)=mean_diff;
        var_of_sum_group(1,tmp)=var_of_sum;
        posterior_probability_group(1,tmp)=posterior_probability;
            
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/GSR/Smith/Full_model/Lateralization_index_A_comp_GSR_DMN.mat']);
        mean_diff_group(2,tmp)=mean_diff;
        var_of_sum_group(2,tmp)=var_of_sum;
        posterior_probability_group(2,tmp)=posterior_probability;
        
        
    end
end

save([Work_dir '/Results_paper_variability/DCM/GSR/Smith/Full_model/Subject_asymmetry.mat'],'mean_diff_group','var_of_sum_group','posterior_probability_group');

%%%%%%%%%%%%%%%%%%%%%%%
%Reliability: GSR
%%%%%%%%%%%%%%%%%%%%%%%
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These subjects were excluded based on diagnostic checks (see previous codes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
            
            
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/WS_correlation_comp_GSR_DMN.mat']);
        cr_group_GSR(1,tmp)=mean(cr_6);
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/GSR/Smith/Full_model/WS_correlation_comp_GSR_DMN.mat']);
        cr_group_GSR(2,tmp)=mean(cr_6);

    end
end

cr_group_GSR(:,end+1)=mean(cr_group_GSR,2);

save([Work_dir '/Results_paper_variability/DCM/GSR/Smith/Full_model/WS_correlation_GSR.mat'],'cr_group_GSR');

end