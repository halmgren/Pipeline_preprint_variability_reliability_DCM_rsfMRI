function Wrapper_extract_timeseries_paper_variability(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to parallelly extract VOIs for several sessions

for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    
    parfor session=1:length(session_number)
        session_name=session_number(session).name;
        Extract_timeseries_paper_variability(dataset,session_name,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
    end
end

end