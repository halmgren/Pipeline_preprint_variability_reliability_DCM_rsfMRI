function PEB_reduced_session_level_paper_variability(SPM_dir,Work_dir)

for number_dataset=1:4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract information about datasets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%02d',subject)]);
        session_number(1:2)=[];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute hemispheric dominance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        name_ROI_def='Smith';
        procedure='BMR';
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
        
        for network_number=1:length(ntwrk_name)
            disp(ntwrk_name{network_number});
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
                
                c=NaN(1,length(DCM));
                v=NaN(1,length(DCM));
                PP=NaN(1,length(DCM));
                for session=1:length(DCM)
                    
                    %DMN
                    if strcmp(ntwrk_name{network_number},'DMN')
                        T = 0;
                        C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(spm_vec(DCM{session}.Ep))-16)]'/3;
                        c(session) = C'*spm_vec(DCM{session}.Ep);
                        v(session) = C'*DCM{session}.Cp*C;
                        PP(session)   = 1-spm_Ncdf(T,c(session),v(session));
                    end
                    
                    
                end
                
                mean_diff=c;
                var_of_sum=v;
                posterior_probability=PP;
                
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_individ_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                clear mean_diff var_of_sum c v C PP T posterior_probability;
                
                clear PEB DCM;
            catch
            end
        end
        
        %%%%%%%%%%%%%%%%%%%
        %Cross-correlation
        %%%%%%%%%%%%%%%%%%%
        for network_number=1:length(ntwrk_name)
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
                
                clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                
                for session=1:length(DCM)-1
                    cr_1(session)=corr(DCM{session}.Ep.A(:),DCM{session+1}.Ep.A(:));
                end
                
                
                for session=1:length(DCM)-2
                    cr_2(session)=corr(DCM{session}.Ep.A(:),DCM{session+2}.Ep.A(:));
                end
                
                
                for session=1:length(DCM)-3
                    cr_3(session)=corr(DCM{session}.Ep.A(:),DCM{session+3}.Ep.A(:));
                end
                
                
                for session=1:length(DCM)-4
                    cr_4(session)=corr(DCM{session}.Ep.A(:),DCM{session+4}.Ep.A(:));
                end
                
                
                for session=1:length(DCM)-5
                    cr_5(session)=corr(DCM{session}.Ep.A(:),DCM{session+5}.Ep.A(:));
                end
                
                B=nchoosek(1:length(DCM),2);
                
                for combi=1:length(B)
                    cr_6(combi)=corr(DCM{B(combi,1)}.Ep.A(:),DCM{B(combi,2)}.Ep.A(:));
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/WS_correlation_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6');
                clear DCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            catch
                continue;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Delat cross-correlation
        %%%%%%%%%%%%%%%%%%%%%%%%%
        for network_number=1:length(ntwrk_name)
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_individ_' ntwrk_name{network_number} '.mat']);
                
                clear cr_1 cr_2 cr_3 cr_4 cr_5 cr_6;
                
                
                DCM_delat=DCM;
                
                for session=1:length(DCM)
                    if mean_diff(session)>0
                        DCM_delat{session}.Ep.A([2 4],:)=DCM_delat{session}.Ep.A([4 2],:);
                        DCM_delat{session}.Ep.A(:,[2 4])=DCM_delat{session}.Ep.A(:,[4 2]);
                    end
                end
                
                for session=1:length(DCM_delat)-1
                    cr_1(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+1}.Ep.A(:));
                end
                
                
                for session=1:length(DCM_delat)-2
                    cr_2(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+2}.Ep.A(:));
                end
                
                
                for session=1:length(DCM_delat)-3
                    cr_3(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+3}.Ep.A(:));
                end
                
                
                for session=1:length(DCM_delat)-4
                    cr_4(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+4}.Ep.A(:));
                end
                
                
                for session=1:length(DCM_delat)-5
                    cr_5(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+5}.Ep.A(:));
                end
                
                
                for session=1:length(DCM_delat)-5
                    cr_5(session)=corr(DCM_delat{session}.Ep.A(:),DCM_delat{session+5}.Ep.A(:));
                end
                
                B=nchoosek(1:length(DCM_delat),2);
                
                for combi=1:length(B)
                    cr_6(combi)=corr(DCM_delat{B(combi,1)}.Ep.A(:),DCM_delat{B(combi,2)}.Ep.A(:));
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/WS_correlation_delat_' ntwrk_name{network_number} '.mat'],'cr_1','cr_2','cr_3','cr_4','cr_5','cr_6')
                clear DCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            catch
                continue;
            end
        end
        
        
    end
end
end