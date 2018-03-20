function PEB_subject_mean_paper_variability(dataset,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

disp(dataset);
disp(num2str(subject));
disp(procedure);

%%%%%%%%%%%%%%%%
%Basic Analyses
%%%%%%%%%%%%%%%%
if strcmp(procedure,'Basic')
    
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
    
    for network_number=1:length(ntwrk_name)
        disp(ntwrk_name{network_number});
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
        
        flag=zeros(1,length(GCM));
        for diagn=1:length(GCM)
            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                flag(diagn)=1;
            end
        end
        
        GCM(find(flag==1))=[];
        
        %PEB settings
        M = struct();
        M.alpha = 1;
        M.beta  = 16;
        M.hE    = 0;
        M.hC    = 1/16;
        M.Q     = 'all';
        M.X=ones(length(GCM),1);
        
        [PEB DCM]=spm_dcm_peb(GCM,M,'A');
        
        save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
        
        clear GCM PEB DCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%
%GSR-Basic comparison
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
    
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
    
    for network_number=1:length(ntwrk_name)
        disp(ntwrk_name{network_number});
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
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
        
        %PEB settings
        M = struct();
        M.alpha = 1;
        M.beta  = 16;
        M.hE    = 0;
        M.hC    = 1/16;
        M.Q     = 'all';
        M.X=ones(length(GCM),1);
        
        
        [PEB DCM]=spm_dcm_peb(GCM,M,'A');
        
        save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
        
        clear GCM PEB DCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROI sizes: Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(procedure,'ROI_Size')
    
    for size_ROI=1:4
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
        
        for network_number=1:length(ntwrk_name)
            
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
            
            %PEB settings
            M = struct();
            M.alpha = 1;
            M.beta  = 16;
            M.hE    = 0;
            M.hC    = 1/16;
            M.Q     = 'all';
            M.X=ones(length(GCM),1);
            
            [PEB DCM]=spm_dcm_peb(GCM,M,'A');
            
            save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
            
            clear GCM PEB DCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute hemispheric dominance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(procedure,'Basic')
    
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
    
    for network_number=1:length(ntwrk_name)
        disp(ntwrk_name{network_number});
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
        
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
        
        save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
        clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute hemispheric dominance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
    
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
    
    for network_number=1:length(ntwrk_name)
        disp(ntwrk_name{network_number});
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
        
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
        
        save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_A_comp_GSR_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
        clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
    end
end
 
if strcmp(procedure,'ROI_Size')
    
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
    
    for network_number=1:length(ntwrk_name)
        for size_ROI=1:4
            disp(ntwrk_name{network_number});
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
            
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
            
            save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/Lateralization_index_A_comp_ROIsize_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
            clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
        end
    end
end

end