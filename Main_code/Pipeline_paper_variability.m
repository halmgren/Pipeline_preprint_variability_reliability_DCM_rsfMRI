%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code for the paper 'Variability and reliabilty of effective connectivity within the core decault mode network'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOFTWARE VERSIONS
%
%OS: ubuntu 16.04 LTS; Intel® Core™ i7-6800K CPU @ 3.40GHz × 12
%MATLAB: 8.6.0 (R2015b)
%SPM: SPM12 (revision 6906)
%DCM: DCM12 (revision 6801)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATASETS
%
%'MyConnectome',            Downloaded from https://openfmri.org/dataset/ds000031/
%'The Midnight Scan Club'   Downloaded from https://openfmri.org/dataset/ds000224/
%'Kirby'                    Downloaded from https://www.nitrc.org/projects/kirbyweekly
%'Day2day'                  For availability, see Filevich et al. (2017); section 'Availability of data and materials'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATASET PREPARATION
%
%Because different datasets were structured in different formats (which might change in future 
%releases), they should all be (re)structured according to the BIDS format (http://bids.neuroimaging.io/).
%In addition, 4D resting state files should be renamed as 'functional.nii' and T1-weighted images should be 
%renamed as 'structural.nii'. In future releases we will try to encorporate some of these (re)structuring steps in
%the code.
%
%The datasets should thus be structured as follows:
%
%Work_dir
%   my_dataset1
%       sub-01          (for subject two digits are required, e.g., sub-05 or sub-10)
%           ses-001     (for session three digits are required, e.g., ses-005 or ses-015)
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ses-002
%               func
%                   functional.nii
%           ...
%
%       sub-02
%           ses-001
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ...
%       ...
%
%   my_dataset2
%       sub-01
%           ses-001
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ...
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!
%THE NAMES OF THE DATASET-FOLDERS ARE FIXED:
%'Myconnectome':            change 'my_dataset' to 'DatasetPoldrack' 
%'The Midnight Scan Club':  change 'my_dataset' to 'DatasetGordon' 
%'Kirby':                   change 'my_dataset' to 'DatasetKirby' 
%'Day2day':                 change 'my_dataset' to 'DatasetKuehn'
%!!!!!!!!!!!!!!!!!!!!!!!!!
%
%Anatomical (T1) images that were used, were stored in the 'anat' folder of the FIRST session that was analyzed
%
%               - 'MyConnectome':         T1 image of session 12 was used and should be stored in 'anat' folder of session 13 (first fMRI session analyzed)
%               - 'Midnight Scan Club':   First T1 image acquired (for each subject) was used and should be stored in 'anat' folder of session 1
%               - 'Day2day' and 'Kirby':  T1 image was acquired on same day as first fMRI scan
%
%All neuroimaging data should be in 4D nifti format AND should be unpacked/extracted
%
%The number of iterations in spm_dcm_peb.m was increased to 128 (changed 'for n = 1:64' to 'for n = 1:128')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURES
%
%For some figures we made use of the 'gramm'-toolbox (release 2.23), downloadable from https://github.com/piermorel/gramm
%
%For other figures we made use of the 'igraph' package (version 1.1.2) AND the 'R.matlab' package (version 3.6.1) in R (R version 3.2.3) 
%   -> make sure you also change the 'Work_dir' file in the R-script to the same as in this script
%
%We also used the robust statistical toolbox (https://github.com/CPernet/Robust_Statistical_Toolbox) AND the scatterHistDiff.m code (https://github.com/anne-urai/Tools/blob/master/plotting/scatterHistDiff.m)
%   -> For both, the adapted versions are included in the present code-package (no need to download separately)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We encountered two duplicates in the 'Kirby' dataset. These are deleted from the present release, so if 
%you have downloaded the dataset before 2018-01-02 you should redownload it before using this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For comments and feedback, please contact Hannes_Almgren@ugent.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify SPM directory and working directory for datasets (do NOT include '/' at the end)
SPM_dir='/home/hannes/Desktop/spm12';
Work_dir='/media/hannes/Almgren_Disk3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%From here on everything is automatic
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

Move_artefact_and_pilot_scans(Work_dir)

Scanning_parameters_paper_variability(Work_dir) %Create files with information regarding scan parameters (TR etc)

Move_MRI_files_paper_variability(Work_dir) %Move functional images to separate folder

Remove_first_scans_paper_variability(Work_dir) %Discard first five scans of each session

for number_dataset=1:4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Information about datasets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009
    all_procedure_names={'Basic','ROI_Size','GSR'};
    
    for number_ROI_def=1:length(all_ROI_defs)   %the code can be adapted to include different ROI definitions (e.g., different coordinates) => if changed, should also be changed in many other codes (e.g., PEB_group_mean_paper_variability)
        
        name_ROI_def=all_ROI_defs{number_ROI_def};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Give regions name and coördinates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %DMN was studied for this paper
        [ROI_list]=Define_ROIs_paper_variability(name_ROI_def);
        
        for number_procedure=1:length(all_procedure_names)
            
            procedure=all_procedure_names{number_procedure};
            
            %%%%%%%%%%%%%%%%
            %Thresholds QC
            %%%%%%%%%%%%%%%%
            [threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_threshold_VOIs]=Define_QC_tresholds_paper_variability(procedure);
            
            %%%%%%%%%%%%%%%%%
            %%Preprocessing
            %%%%%%%%%%%%%%%%%
            %structural
            Preprocess_structural_paper_variability(dataset,number_subject,procedure,SPM_dir,Work_dir);
            
            %Functional
            Wrapper_preprocess_functional_paper_variability(dataset,number_subject,slice_time_seconds,procedure,SPM_dir,Work_dir);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Time-series extraction
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            %Extract regressors for CSV, WM, and discrete cosine set
            Wrapper_extract_regressors_paper_variability(dataset,number_subject,procedure,SPM_dir,Work_dir);
            
            %Extract timeseries for all regions
            Wrapper_extract_timeseries_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Specify and estimate DCMs: Full model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Separate estimation
            Wrapper_specify_estimate_DCM_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)
            
            %%Save diagnostics such as number of estimable parameters, explained variance, max connection strength, excessive motion, alpha level ROI extraction, etc...
            Diagnostics_paper_variability(dataset,number_subject,ROI_list,name_ROI_def,procedure,SPM_dir,Work_dir)
            
            %Compute and save whether session reach specifie threshold
            Above_threshold_sessions_paper_variability(dataset,number_subject,ROI_list,name_ROI_def,threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_threshold_VOIs,procedure,SPM_dir,Work_dir)
            
        end

    end
end


for number_dataset=1:4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Information about datasets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009
    all_procedure_names={'Basic','ROI_Size','GSR'};
    
    for number_ROI_def=1:length(all_ROI_defs)   %the code can be adapted to include different ROI definitions (e.g., different coordinates) => if changed, should also be changed in many other codes (e.g., PEB_group_mean_paper_variability)
        
        name_ROI_def=all_ROI_defs{number_ROI_def};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Give regions name and coördinates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %DMN was studied for this paper
        [ROI_list]=Define_ROIs_paper_variability(name_ROI_def);
        
        for number_procedure=1:length(all_procedure_names)
            
            procedure=all_procedure_names{number_procedure};
            
            %%%%%%%%%%%%%%%%
            %Thresholds QC
            %%%%%%%%%%%%%%%%
            [threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_threshold_VOIs]=Define_QC_tresholds_paper_variability(procedure);
            
            %%%%%%%%%%
            %PEB Mean
            %%%%%%%%%%%
            Wrapper_PEB_subject_mean_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)
            
            Compute_WS_correlation_paper_variability(dataset,number_subject,procedure,name_ROI_def,ROI_list,SPM_dir,Work_dir)
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Summary across subjects
%%%%%%%%%%%%%%%%%%%%%%%%%
PEB_group_mean_paper_variability(SPM_dir,Work_dir);

Compute_BS_correlation_paper_variability(SPM_dir,Work_dir);

Compute_BS_sign_stability_paper_variability(SPM_dir,Work_dir);

PEB_reduced_session_level_paper_variability(SPM_dir,Work_dir);

PEB_reduced_subject_level_paper_variability(SPM_dir,Work_dir);

PCA_BS_paper_variability(SPM_dir,Work_dir);

PCA_WS_paper_variability(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compare parameter structures between and within procedures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Compute_between_method_correlation_paper_variability(SPM_dir,Work_dir);

%%%%%%%%%
%Figures
%%%%%%%%%
Figures_paper_variability(SPM_dir,Work_dir);

%%%%%%%%%
%Results
%%%%%%%%%
Results_paper_variability(SPM_dir,Work_dir);
