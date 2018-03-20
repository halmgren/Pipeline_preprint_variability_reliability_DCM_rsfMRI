function Compute_WS_correlation_paper_variability(dataset,number_subject,procedure,name_ROI_def,ROI_list,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only cr_6 (all permutations) were analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute within-subject correlation for original matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subject=1:number_subject
    
    if strcmp(procedure,'Basic')
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model'])
        
        for network_number=1:length(ntwrk_name)
            
            disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
            try
                load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                
                flag=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag(diagn)=1;
                    end
                end
                
                GCM(find(flag==1))=[];
                
                clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                
                for session=1:length(GCM)-1
                    cr_1(session)=corr(GCM{session}.Ep.A(:),GCM{session+1}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-2
                    cr_2(session)=corr(GCM{session}.Ep.A(:),GCM{session+2}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-3
                    cr_3(session)=corr(GCM{session}.Ep.A(:),GCM{session+3}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-4
                    cr_4(session)=corr(GCM{session}.Ep.A(:),GCM{session+4}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-5
                    cr_5(session)=corr(GCM{session}.Ep.A(:),GCM{session+5}.Ep.A(:));
                end
                
                B=nchoosek(1:length(GCM),2);
                
                for combi=1:length(B)
                    cr_6(combi)=corr(GCM{B(combi,1)}.Ep.A(:),GCM{B(combi,2)}.Ep.A(:));
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/WS_correlation_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6')
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            catch ME
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_WS_correlation_' ntwrk_name{network_number} '.mat'],'ME');
                continue
            end
        end
        
    end
    
    if strcmp(procedure,'GSR')||strcmp(procedure,'Basic')
        
        for network_number=1:length(ntwrk_name)
            try
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                
                flag=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag(diagn)=1;
                    end
                end
                
                if strcmp(procedure,'Basic')
                    procedure_comp='GSR';
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                end
                
                if strcmp(procedure,'GSR')
                    procedure_comp='Basic';
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                end
                
                flag_comp=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag_comp(diagn)=1;
                    end
                end
                
                GCM(find(flag==1|flag_comp==1))=[];
                
                clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                
                for session=1:length(GCM)-1
                    cr_1(session)=corr(GCM{session}.Ep.A(:),GCM{session+1}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-2
                    cr_2(session)=corr(GCM{session}.Ep.A(:),GCM{session+2}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-3
                    cr_3(session)=corr(GCM{session}.Ep.A(:),GCM{session+3}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-4
                    cr_4(session)=corr(GCM{session}.Ep.A(:),GCM{session+4}.Ep.A(:));
                end
                
                
                for session=1:length(GCM)-5
                    cr_5(session)=corr(GCM{session}.Ep.A(:),GCM{session+5}.Ep.A(:));
                end
                
                B=nchoosek(1:length(GCM),2);
                
                for combi=1:length(B)
                    cr_6(combi)=corr(GCM{B(combi,1)}.Ep.A(:),GCM{B(combi,2)}.Ep.A(:));
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/WS_correlation_comp_GSR_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6');
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            catch ME
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_WS_correlation_GSR_' ntwrk_name{network_number} '.mat'],'ME');
                continue
            end
        end
        
    end
    
    if strcmp(procedure,'ROI_Size')
        for size_ROI=1:4
            
            for network_number=1:length(ntwrk_name)
                try
                    cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                    
                    flag=zeros(1,length(GCM));
                    for diagn=1:length(GCM)
                        if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                            flag(diagn)=1;
                        end
                    end
                    
                    R=1:4;
                    R(size_ROI)=[];
                    
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(1)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                    
                    flag_comp1=zeros(1,length(GCM));
                    for diagn=1:length(GCM)
                        if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                            flag_comp1(diagn)=1;
                        end
                    end
                    
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(2)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                    
                    flag_comp2=zeros(1,length(GCM));
                    for diagn=1:length(GCM)
                        if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                            flag_comp2(diagn)=1;
                        end
                    end
                    
                    load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(3)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                    
                    flag_comp3=zeros(1,length(GCM));
                    for diagn=1:length(GCM)
                        if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                            flag_comp3(diagn)=1;
                        end
                    end
                    
                    
                    GCM(find(flag==1|flag_comp1==1|flag_comp2==1|flag_comp3==1))=[];
                    
                    clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                    
                    for session=1:length(GCM)-1
                        cr_1(session)=corr(GCM{session}.Ep.A(:),GCM{session+1}.Ep.A(:));
                    end
                    
                    
                    for session=1:length(GCM)-2
                        cr_2(session)=corr(GCM{session}.Ep.A(:),GCM{session+2}.Ep.A(:));
                    end
                    
                    
                    for session=1:length(GCM)-3
                        cr_3(session)=corr(GCM{session}.Ep.A(:),GCM{session+3}.Ep.A(:));
                    end
                    
                    
                    for session=1:length(GCM)-4
                        cr_4(session)=corr(GCM{session}.Ep.A(:),GCM{session+4}.Ep.A(:));
                    end
                    
                    
                    for session=1:length(GCM)-5
                        cr_5(session)=corr(GCM{session}.Ep.A(:),GCM{session+5}.Ep.A(:));
                    end
                    
                    B=nchoosek(1:length(GCM),2);
                    
                    for combi=1:length(B)
                        cr_6(combi)=corr(GCM{B(combi,1)}.Ep.A(:),GCM{B(combi,2)}.Ep.A(:));
                    end
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/WS_correlation_comp_ROIsize_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6');
                    clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
                catch ME
                    save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/ERROR_WS_correlation_' ntwrk_name{network_number} '.mat'],'ME');
                    continue
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reorder matrices such that dominant hemispheres are at same side
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(procedure,'Basic')
        
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model'])
        
        for network_number=1:length(ntwrk_name)
            disp(['Dataset: ' dataset '; subject: ' num2str(subject) '; network: ' ntwrk_name{network_number}]);
            try
                load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateral_index_individ_' ntwrk_name{network_number} '.mat']);
                
                flag=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag(diagn)=1;
                    end
                end
                
                GCM_delat=GCM;
                
                clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                
                for session=1:length(GCM)
                    if mean_diff(session)>0
                        GCM_delat{session}.Ep.A([2 4],:)=GCM_delat{session}.Ep.A([4 2],:);
                        GCM_delat{session}.Ep.A(:,[2 4])=GCM_delat{session}.Ep.A(:,[4 2]);
                    end
                end
                
                GCM_delat(find(flag==1))=[];
                
                for session=1:length(GCM_delat)-1
                    cr_1(session)=corr(GCM_delat{session}.Ep.A(:),GCM_delat{session+1}.Ep.A(:));
                end
                
                
                for session=1:length(GCM_delat)-2
                    cr_2(session)=corr(GCM_delat{session}.Ep.A(:),GCM_delat{session+2}.Ep.A(:));
                end
                
                
                for session=1:length(GCM_delat)-3
                    cr_3(session)=corr(GCM_delat{session}.Ep.A(:),GCM_delat{session+3}.Ep.A(:));
                end
                
                
                for session=1:length(GCM_delat)-4
                    cr_4(session)=corr(GCM_delat{session}.Ep.A(:),GCM_delat{session+4}.Ep.A(:));
                end
                
                
                for session=1:length(GCM_delat)-5
                    cr_5(session)=corr(GCM_delat{session}.Ep.A(:),GCM_delat{session+5}.Ep.A(:));
                end
                
                B=nchoosek(1:length(GCM_delat),2);
                
                for combi=1:length(B)
                    cr_6(combi)=corr(GCM_delat{B(combi,1)}.Ep.A(:),GCM_delat{B(combi,2)}.Ep.A(:));
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/WS_correlation_delat_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6')
                clear GCM GCM_delat Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            catch ME
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_WS_correlation_delat_' ntwrk_name{network_number} '.mat'],'ME');
                continue
            end
        end
    end
    
end
end