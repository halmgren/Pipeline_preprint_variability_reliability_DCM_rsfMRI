function Wrapper_extract_regressors_paper_variability(dataset,number_subject,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to  extract regressors for several sessions parallelly

for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    
    parfor session=1:length(session_number)
        session_name=session_number(session).name;
        Extract_regressors_paper_variability(dataset,session_name,subject,procedure,SPM_dir,Work_dir);
    end
end

end