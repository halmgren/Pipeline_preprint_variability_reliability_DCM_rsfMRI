function Specify_estimate_DCM_paper_variability(dataset,session_name,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%
%specify spDCM file
%%%%%%%%%%%%%%%%%%%%%
disp('Specify and estimate full model DCM')
disp(['Dataset: ' dataset '; subject: ' num2str(subject) '; session: ' session_name]);

if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
    cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);
    
    tmp=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            continue
        else
            tmp=tmp+1;
            ntwrk_name{tmp}=spm_select('FPList',pwd,['^VOI_' ntwrk '.*\.mat$']);
            
        end
    end
    
    cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/scan_info/']);
    load('General_information.mat');
    
    cd('../fourD_files/Basic');
    [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);    %!!!only used to compute number of scans!!!
    length_scan=size(nifti_images,1);
    
    for network_number=1:length(ntwrk_name)
        clear DCM;
        number_regions=size(ntwrk_name{network_number},1);
        
        DCM.a=ones(number_regions,number_regions);
        DCM.b=zeros(number_regions,number_regions);
        DCM.c=zeros(number_regions,1);
        DCM.d=zeros(number_regions,number_regions,0);
        
        DCM.U.u=zeros(length_scan,1);
        DCM.U.name={'null'};
        
        xY=[];
        for region_number=1:number_regions
            if isspace(ntwrk_name{network_number}(region_number,end))==1
                K=load(ntwrk_name{network_number}(region_number,1:end-1),'xY');
            else
                K=load(ntwrk_name{network_number}(region_number,1:end),'xY');
            end
            xY = spm_cat_struct(xY,K.xY);
        end
        
        DCM.xY=xY;
        
        DCM.Y.dt=TR;
        DCM.Y.X0=xY(1).X0;
        DCM.Y.Q=spm_Ce(ones(1,number_regions)*length_scan);
        
        for  region_number = 1:number_regions
            DCM.Y.y(:,region_number)  = xY(region_number).u;
            DCM.Y.name{region_number} = xY(region_number).name;
        end
        
        DCM.v=length_scan;
        DCM.n=number_regions;
        
        DCM.TE=TE;
        DCM.delays=ones(number_regions,1)*TR/2;
        
        DCM.options.nonlinear=0;
        DCM.options.two_state=0;
        DCM.options.stochastic=0;
        DCM.options.centre=0;
        DCM.options.induced=1;
        disp(num2str(number_regions))
        GCM_full{network_number}=DCM;
    end
    
    GCM_full=GCM_full'; %column cell array
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Estimate DCM files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
    
    cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
    
    GCM=spm_dcm_fit(GCM_full);
    save('GCM_full_estim.mat','GCM');
    GCM=GCM_full;
    save('GCM_full_raw.mat','GCM');
    delete([sprintf('DCM_%s',date) '.mat']);
    clear GCM;
    
elseif strcmp(procedure,'ROI_Size')
    for size_ROI=1:4
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def '/' num2str(4*size_ROI) '/']);
        tmp=0;
        
        for VOI_number=1:size(ROI_list,1)
            ntwrk=ROI_list{VOI_number,1}(1:3);
            
            if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
                continue
            else
                tmp=tmp+1;
                ntwrk_name{tmp}=spm_select('FPList',pwd,['^VOI_' ntwrk '.*\.mat$']);
            end
        end
        
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/scan_info/']);
        load('General_information.mat');
        
        cd('../fourD_files/Basic');
        [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);
        length_scan=size(nifti_images,1);
        
        for network_number=1:length(ntwrk_name)
            clear DCM;
            number_regions=size(ntwrk_name{network_number},1);
            
            DCM.a=ones(number_regions,number_regions);
            DCM.b=zeros(number_regions,number_regions);
            DCM.c=zeros(number_regions,1);
            DCM.d=zeros(number_regions,number_regions,0);
            
            DCM.U.u=zeros(length_scan,1);
            DCM.U.name={'null'};
            
            xY=[];
            for region_number=1:number_regions
                if isspace(ntwrk_name{network_number}(region_number,end))==1
                    K=load(ntwrk_name{network_number}(region_number,1:end-1),'xY');
                else
                    K=load(ntwrk_name{network_number}(region_number,1:end),'xY');
                end
                xY = spm_cat_struct(xY,K.xY);
            end
            
            DCM.xY=xY;
            
            DCM.Y.dt=TR;
            DCM.Y.X0=xY(1).X0;
            DCM.Y.Q=spm_Ce(ones(1,number_regions)*length_scan);
            
            for  region_number = 1:number_regions
                DCM.Y.y(:,region_number)  = xY(region_number).u;
                DCM.Y.name{region_number} = xY(region_number).name;
            end
            
            DCM.v=length_scan;
            DCM.n=number_regions;
            
            DCM.TE=TE;
            DCM.delays=ones(number_regions,1)*TR/2;
            
            DCM.options.nonlinear=0;
            DCM.options.two_state=0;
            DCM.options.stochastic=0;
            DCM.options.centre=0;
            DCM.options.induced=1;
            disp(num2str(number_regions))
            GCM_full{network_number}=DCM;
        end
        
        GCM_full=GCM_full'; %column cell array
        
        mkdir([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
        
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)]);
        
        GCM=spm_dcm_fit(GCM_full);
        save('GCM_full_estim.mat','GCM');
        GCM=GCM_full;
        save('GCM_full_raw.mat','GCM');
        delete([sprintf('DCM_%s',date) '.mat']);
        clear GCM;
        
    end
    
end

end