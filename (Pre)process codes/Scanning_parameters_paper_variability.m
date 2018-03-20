function Scanning_parameters_paper_variability(Work_dir)

    %%%%%%%%%%%%%%%%
    %Day2day
    %%%%%%%%%%%%%%%%
    image_dim=[72,72,36];
    interslice_gap=0.6;
    N_scans=150;
    TE=0.03;
    TR=2;
    vox_size=[3,3,3];
    
    slice_times=1.0e+03 * [1.0025, 0, 1.0600, 0.0550, 1.1150, 0.1100, 1.1700, 0.1675, 1.2275, 0.2225, 1.2825, 0.2775, 1.3375, 0.3350, 1.3925, 0.3900, 1.4500, 0.4450, 1.5050, 0.5025, 1.5600, 0.5575, 1.6175, 0.6125, 1.6725, 0.6700, 1.7275, 0.7250, 1.7850, 0.7800, 1.8400, 0.8350, 1.8950, 0.8925, 1.9525, 0.9475];
    
    dataset='DatasetKuehn';
    for subject=1:8
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('scan_info');
            save('scan_info/General_information.mat','image_dim','interslice_gap','N_scans','TE','TR','vox_size');
            save('scan_info/slice_timing.mat','slice_times');
            cd ../..;
        end
        
    end
    
    clear image_dim interlice_gap N_scans TE TR vox_size slice_times;
    
    %%%%%%%%%%%%%%%%
    %MyConnectome
    %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%
    %Session 13 till 26
    %%%%%%%%%%%%%%%%%%%%
    
    
    image_dim=[96 96 68];
    interslice_gap=0.4;
    N_scans=518;
    TE=0.03;
    TR=1.16;
    vox_size=[2.4,2.4,2];
    
    slice_times=[0,615,68,682,138,750,205,820,273,888,342,955,410,1023,478,1092,545,0,615,68,682,138,750,205,820,273,888,342,955,410,1023,478,1092,545,0,615,68,682,138,750,205,820,273,888,342,955,410,1023,478,1092,545,0,615,68,682,138,750,205,820,273,888,342,955,410,1023,478,1092,545];
    
    dataset='DatasetPoldrack';
    cd([Work_dir '/' dataset '/sub-01']);
    session_number=dir;
    session_number(1:2)=[];
    
    for session=1:14
        session_name=session_number(session).name;
        cd([session_name '/func']);
        mkdir('scan_info');
        save('scan_info/General_information.mat','image_dim','interslice_gap','N_scans','TE','TR','vox_size');
        save('scan_info/slice_timing.mat','slice_times');
        cd ../..;
    end
    
    clear image_dim interlice_gap N_scans TE TR vox_size slice_times;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Session 27 and further
    %%%%%%%%%%%%%%%%%%%%%%%%
     
    image_dim=[96 96 64];
    interslice_gap=0.4;
    N_scans=518;
    TE=0.03;
    TR=1.16;
    vox_size=[2.4,2.4,2];
    
    slice_times=[0,498,995,355,853,213,710,70,568,1068,425,925,285,782,143,640,0,498,995,355,853,213,710,70,568,1068,425,925,285,782,143,640,0,498,995,355,853,213,710,70,568,1068,425,925,285,782,143,640,0,498,995,355,853,213,710,70,568,1068,425,925,285,782,143,640];
    
    dataset='DatasetPoldrack';
    cd([Work_dir '/' dataset '/sub-01']);
    session_number=dir;
    session_number(1:2)=[];
    
    for session=15:length(session_number)
        session_name=session_number(session).name;
        cd([session_name '/func']);
        mkdir('scan_info');
        save('scan_info/General_information.mat','image_dim','interslice_gap','N_scans','TE','TR','vox_size');
        save('scan_info/slice_timing.mat','slice_times');
        cd ../..;
    end
    
    clear image_dim interlice_gap N_scans TE TR vox_size slice_times;
    
    %%%%%%%%%%%%%%%%
    %Kirby
    %%%%%%%%%%%%%%%%
    
    image_dim=[80 80 37];
    interslice_gap=1;
    N_scans=200;
    TE=0.03;
    TR=2;
    vox_size=[3 3 3];
    
    slice_times=[1:37]; %ascending serial
    
    dataset='DatasetKirby';
    cd([Work_dir '/' dataset '/sub-01']);
    session_number=dir;
    session_number(1:2)=[];
    
    for session=1:length(session_number)
        session_name=session_number(session).name;
        cd([session_name '/func']);
        mkdir('scan_info');
        save('scan_info/General_information.mat','image_dim','interslice_gap','N_scans','TE','TR','vox_size');
        save('scan_info/slice_timing.mat','slice_times');
        cd ../..;
    end
    
    clear image_dim interlice_gap N_scans TE TR vox_size slice_times;
    
    %%%%%%%%%%%%%%%%
    %MSC
    %%%%%%%%%%%%%%%%
    
    image_dim=[64 64 36];
    interslice_gap='unknown'; %interslice gap was not used in any analyses
    N_scans=818;
    TE=0.027;
    TR=2.2;
    vox_size=[4,4,4];
    
    slice_times=[1.105, 0, 1.165, 0.06,	1.2275,	0.1225,	1.2875,	0.185, 1.35, 0.245,	1.41, 0.3075, 1.4725, 0.3675, 1.5325, 0.43,	1.595, 0.49, 1.655,	0.5525,	1.7175,	0.6125,	1.78, 0.675, 1.84, 0.735, 1.9025, 0.7975, 1.9625, 0.8575, 2.025, 0.92, 2.085, 0.9825, 2.1475, 1.0425];

    dataset='DatasetGordon';
    
    for subject=1:10
        cd([Work_dir '/' dataset '/sub-' sprintf('%02d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('scan_info');
            save('scan_info/General_information.mat','image_dim','interslice_gap','N_scans','TE','TR','vox_size');
            save('scan_info/slice_timing.mat','slice_times');
            cd ../..;
        end
        
    end
   
    clear image_dim interlice_gap N_scans TE TR vox_size slice_times;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For slice timing, the slice times of the first functional image of the first subject of the respective dataset were included in the present script (except for MyConnectome, for which it was changed
    %from session 15 on). For the original analyses, session-specific slice times were used, however, these differed neglibly from the reference (first) functional image (<5ms)
    %slice timing. Although we are sure it won't make a difference, we hope to include the (very precise) session-specific in this script in the future.
    %This is only applicable to the 'myconnectome' and 'day2day' dataset, since they provided session-specific slice timings.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end