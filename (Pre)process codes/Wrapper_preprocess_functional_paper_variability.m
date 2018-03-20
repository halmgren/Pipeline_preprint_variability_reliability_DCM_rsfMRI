function Wrapper_preprocess_functional_paper_variability(dataset,number_subject,slice_time_seconds,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to parallelly preprocess several sessions

for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    session_one=session_number(1).name;
    
    parfor session=1:length(session_number)
        session_name=session_number(session).name;
        Preprocess_functional_paper_variability(dataset,slice_time_seconds,session_name,subject,session_one,procedure,SPM_dir,Work_dir);
    end
end

end