function Move_MRI_files_paper_variability(Work_dir)

    %%%%%%%%%%%%%%%%%
    %Day2day-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetKuehn';
    for subject=1:8
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end

    %%%%%%%%%%%%%%%%%
    %myconnectome-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetPoldrack';
    for subject=1
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end

    %%%%%%%%%%%%%%%%%
    %Kirby-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetKirby';
    for subject=1
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end
 
    %%%%%%%%%%%%%%%%%
    %MSC-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetGordon';
    for subject=1:10
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end
end