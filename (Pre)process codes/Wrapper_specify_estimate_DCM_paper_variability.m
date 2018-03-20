function Wrapper_specify_estimate_DCM_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to parallelly specify and estimate several sessions
for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    
    %invert all sessions seperately
    parfor session=1:length(session_number)
        session_name=session_number(session).name;
        Specify_estimate_DCM_paper_variability(dataset,session_name,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
    end
    
    %%put all sessions for all networks in a separate file
    %Detect separate networks
    tmp=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            continue
            
        else
            tmp=tmp+1;
            ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
        mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model']);
        
        for network=1:length(ntwrk_name)
            
            for session=1:length(session_number)
                session_name=session_number(session).name;
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                load('GCM_full_estim.mat');
                GCM_full{session}=GCM{network};
            end
            
            %Make network-specific files
            GCM=GCM_full';
            cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
            save(['GCM_' ntwrk_name{network} '_full_estim.mat'],'GCM');
            
            
            clear GCM GCM_full;
        end
    end
    
    if strcmp(procedure,'ROI_Size')
        for size_ROI=1:4
            
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
            
            for network=1:length(ntwrk_name)
                
                for session=1:length(session_number)
                    session_name=session_number(session).name;
                    cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
                    load('GCM_full_estim.mat');
                    GCM_full{session}=GCM{network};
                end
                
                %Make network-specific files
                GCM=GCM_full';
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
                save(['GCM_' ntwrk_name{network} '_full_estim.mat'],'GCM');
                
                clear GCM GCM_full;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Lateralization (adapted from spm_dcm_review.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
        mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/'])
        for network_number=1:length(ntwrk_name)
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
            
            c=NaN(1,length(GCM));
            v=NaN(1,length(GCM));
            PP=NaN(1,length(GCM));
            for session=1:length(GCM)
                
                
                %DMN
                if strcmp(ntwrk_name{network_number},'DMN')
                    T = 0;
                    C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(spm_vec(GCM{session}.Ep))-16)]'/3;
                    c(session) = C'*spm_vec(GCM{session}.Ep);   
                    v(session) = C'*GCM{session}.Cp*C;
                    PP(session)   = 1-spm_Ncdf(T,c(session),v(session));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %For visualization purposes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ntwrk_name{network_number},'DMN')
                    T = 0;
                    C1 = [0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 zeros(1,length(spm_vec(GCM{session}.Ep))-16)]'/3;
                    C2 = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 zeros(1,length(spm_vec(GCM{session}.Ep))-16)]'/3;
                    
                    c1(session) = C1'*spm_vec(GCM{session}.Ep);
                    c2(session) = C2'*spm_vec(GCM{session}.Ep);
                    
                    v1(session) = C1'*GCM{session}.Cp*C1;
                    v2(session) = C2'*GCM{session}.Cp*C2;
                    
                    PP1(session)   = 1-spm_Ncdf(T,c1(session),v1(session));
                    PP2(session)   = 1-spm_Ncdf(T,c2(session),v2(session));
                end
                
            end
            
            mean_diff=c;
            var_of_sum=v;
            posterior_probability=PP;
            if strcmp(ntwrk_name{network_number},'DMN')
                
                mean1=c1;
                var_of_sum1=v1;
                posterior_probability1=PP1;
                
                mean2=c2;
                var_of_sum2=v2;
                posterior_probability2=PP2;
            end
            
            
            save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateral_index_individ_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For visualization purposes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(ntwrk_name{network_number},'DMN')
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Left_influence_individ_' ntwrk_name{network_number} '.mat'],'mean1','var_of_sum1','posterior_probability1');
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Right_influence_individ_' ntwrk_name{network_number} '.mat'],'mean2','var_of_sum2','posterior_probability2');
                clear mean_diff var_of_sum c v C PP T posterior_probability GCM;
            end
        end
    end
    
    if strcmp(procedure,'ROI_Size')
        for size_ROI=1:4
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
            for network_number=1:length(ntwrk_name)
                load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
                
                c=NaN(1,length(GCM));
                v=NaN(1,length(GCM));
                PP=NaN(1,length(GCM));
                for session=1:length(GCM)
                    
                    %DMN
                    if strcmp(ntwrk_name{network_number},'DMN')
                        T = 0;
                        C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(spm_vec(GCM{session}.Ep))-16)]'/3;
                        c(session) = C'*spm_vec(GCM{session}.Ep);
                        v(session) = C'*GCM{session}.Cp*C;
                        PP(session)   = 1-spm_Ncdf(T,c(session),v(session));
                    end
                    
                end
                
                mean_diff=c;
                var_of_sum=v;
                posterior_probability=PP;
                
                save([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/Lateral_index_individ_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                clear mean_diff var_of_sum c v C PP T posterior_probability GCM;
            end
        end
    end
    
end

end