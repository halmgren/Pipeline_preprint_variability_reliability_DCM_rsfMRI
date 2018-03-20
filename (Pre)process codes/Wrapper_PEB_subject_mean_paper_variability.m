function Wrapper_PEB_subject_mean_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

parfor subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    
    PEB_subject_mean_paper_variability(dataset,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);

end

end