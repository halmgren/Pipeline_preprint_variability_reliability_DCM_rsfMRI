function Remove_first_scans_paper_variability(Work_dir)

    %%%%%%%%%%%%%%%%%
    %Day2day-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetKuehn';
    for subject=1:8
        
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            clear matlabbatch;
            session_name=session_number(session).name;
            cd([session_name '/func/fourD_files']);
            [nifti_images,~]=spm_select('ExtFPList',pwd,'^functional_4D.nii$',inf);
            
            matlabbatch{1}.spm.util.cat.vols = cellstr(nifti_images(6:end,:));
            matlabbatch{1}.spm.util.cat.name = 'functional_disc_4D.nii';
            matlabbatch{1}.spm.util.cat.dtype = 0;
            
            spm_jobman('run',matlabbatch);
            cd ../../..;
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
            clear matlabbatch;
            session_name=session_number(session).name;
            cd([session_name '/func/fourD_files']);
            [nifti_images,~]=spm_select('ExtFPList',pwd,'^functional_4D.nii$',inf);
            
            matlabbatch{1}.spm.util.cat.vols = cellstr(nifti_images(6:end,:));
            matlabbatch{1}.spm.util.cat.name = 'functional_disc_4D.nii';
            matlabbatch{1}.spm.util.cat.dtype = 0;
            
            spm_jobman('run',matlabbatch);
            cd ../../..;
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
            clear matlabbatch;
            session_name=session_number(session).name;
            cd([session_name '/func/fourD_files']);
            [nifti_images,~]=spm_select('ExtFPList',pwd,'^functional_4D.nii$',inf);
            
            matlabbatch{1}.spm.util.cat.vols = cellstr(nifti_images(6:end,:));
            matlabbatch{1}.spm.util.cat.name = 'functional_disc_4D.nii';
            matlabbatch{1}.spm.util.cat.dtype = 0;
            
            spm_jobman('run',matlabbatch);
            cd ../../..;
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
            clear matlabbatch;
            session_name=session_number(session).name;
            cd([session_name '/func/fourD_files']);
            [nifti_images,~]=spm_select('ExtFPList',pwd,'^functional_4D.nii$',inf);
            
            matlabbatch{1}.spm.util.cat.vols = cellstr(nifti_images(6:end,:));
            matlabbatch{1}.spm.util.cat.name = 'functional_disc_4D.nii';
            matlabbatch{1}.spm.util.cat.dtype = 0;
            
            spm_jobman('run',matlabbatch);
            cd ../../..;
        end
        
    end
    
    
end