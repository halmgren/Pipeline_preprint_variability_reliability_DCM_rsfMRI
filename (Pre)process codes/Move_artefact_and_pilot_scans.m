function Move_artefact_and_pilot_scans(Work_dir)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Deleted because of artifacts or related problems
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    movefile([Work_dir '/DatasetKuehn/sub-08/ses-047'],[Work_dir '/DatasetKuehn/sub-08_ses-047_removed']); %Too few scans
    
    movefile([Work_dir '/DatasetKirby/sub-01/ses-012'],[Work_dir '/DatasetKirby/sub-01_ses-012_removed']); %Unsolvable problems in normalization step
    
    movefile([Work_dir '/DatasetPoldrack/sub-01/ses-094'],[Work_dir '/DatasetPoldrack/sub-01_ses-094_removed']); %Low SNR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pilot Sessions 'MyConnectome' study (Not considered in number of rsfMRI sessions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session=1:12
        movefile([Work_dir '/DatasetPoldrack/sub-01/ses-' sprintf('%03d',session)],[Work_dir '/DatasetPoldrack/sub-01_ses-' sprintf('%03d',session) '_removed']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Not considered in number of rsfMRI sessions because (1) not included or (2) not able to open
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    movefile([Work_dir '/DatasetPoldrack/sub-01/ses-052'],[Work_dir '/DatasetPoldrack/sub-01_ses-052_removed']); %No rsfMRI image
    
    movefile([Work_dir '/DatasetPoldrack/sub-01/ses-090'],[Work_dir '/DatasetPoldrack/sub-01_ses-090_removed']); %No rsfMRI image
    
    movefile([Work_dir '/DatasetPoldrack/sub-01/ses-093'],[Work_dir '/DatasetPoldrack/sub-01_ses-093_removed']); %Unable to open and extract rsfMRI
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Two scans sessions were not included in the shared dataset; this code is for possible (future) releases that might include these sessions (Not considered in number of rsfMRI sessions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([Work_dir '/DatasetKuehn/sub-05/ses-043'],'dir')==7 
        movefile([Work_dir '/DatasetKuehn/sub-05/ses-043'],[Work_dir '/DatasetKuehn/sub-05_ses-043_removed']);
    end
    
    if exist([Work_dir '/DatasetKuehn/sub-08/ses-033'],'dir')==7
        movefile([Work_dir '/DatasetKuehn/sub-08/ses-033'],[Work_dir '/DatasetKuehn/sub-08_ses-033_removed']);
        
    end
    
end
