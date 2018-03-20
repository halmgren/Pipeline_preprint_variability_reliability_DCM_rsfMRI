function Preprocess_functional_paper_variability(dataset,slice_time_seconds,session_name,subject,session_one,procedure,SPM_dir,Work_dir)

disp(['subject: ' num2str(subject) ' session: ' session_name]);
disp('Preprocess_functional');
cd([Work_dir '/' dataset '/']);
if strcmp(procedure,'Basic') %Same for all procedure, so only done once
    
    clear matlabbatch;
    
    %extract information regarding TR and slice timing:
    
    cd(['sub-' sprintf('%02d',subject) '/' session_name '/func']);
    cd('scan_info');
    load('General_information.mat');
    load('slice_timing.mat');
    
    % specify path to nifti files
    cd('../fourD_files');
    [nifti_images,~]=spm_select('ExtFPList',pwd,'^functional_disc.*\.nii$',inf);
    
    %Slice time correction information
    matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(nifti_images);
    matlabbatch{1}.spm.temporal.st.nslices = image_dim(3);
    matlabbatch{1}.spm.temporal.st.tr = TR;
    matlabbatch{1}.spm.temporal.st.ta = (slice_time_seconds-1)*-1*(TR-(TR/image_dim(3))); %ignored if so and refslice in ms
    matlabbatch{1}.spm.temporal.st.so = slice_times;
    
    if strcmp(dataset,'DatasetKuehn')
        matlabbatch{1}.spm.temporal.st.refslice = slice_times(36);
    elseif strcmp(dataset,'DatasetPoldrack')
        if image_dim(3)==68
            matlabbatch{1}.spm.temporal.st.refslice = slice_times(17);
        elseif image_dim(3)==64
            matlabbatch{1}.spm.temporal.st.refslice = slice_times(2);
        end
    elseif strcmp(dataset,'DatasetKirby')
        matlabbatch{1}.spm.temporal.st.refslice = slice_times(19);
    elseif strcmp(dataset,'DatasetGordon')
        matlabbatch{1}.spm.temporal.st.refslice = slice_times(36);
    end
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %Spatial realignment
    %%%%%%%%%%%%%%%%%%%%%%%
    
    matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %Coregistration (to first anatomical scan)
    %%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = {[Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_one '/anat/' procedure '/Skullstr_biascor_structural.nii,1']};
    matlabbatch{3}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %Normalize
    %%%%%%%%%%%%%%%%%%%%%%%
    
    matlabbatch{4}.spm.spatial.normalise.write.subj.def = {[Work_dir '/' dataset '/sub-' sprintf('%02d',subject) '/' session_one '/anat/' procedure '/y_structural.nii']};
    matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    tmp=spm_vol(nifti_images(1,:)); %Determine session-specific size for normalisation from functional images (adapted from spm_image.m)
    imat=spm_imatrix(tmp.mat);
    matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = abs(imat(7:9));
    matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %Smooth
    %%%%%%%%%%%%%%%%%%%%%%%
    
    matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{5}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{5}.spm.spatial.smooth.dtype = 0;
    matlabbatch{5}.spm.spatial.smooth.im = 0;
    matlabbatch{5}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run',matlabbatch);
    
    mkdir(procedure);

    movefile('afunctional_disc_4D.mat',procedure);
    movefile('afunctional_disc_4D.nii',procedure);
    movefile('rafunctional_disc_4D.mat',procedure);
    movefile('rafunctional_disc_4D.nii',procedure);
    movefile('wrafunctional_disc_4D.nii',procedure);
    movefile('swrafunctional_disc_4D.nii',procedure);
    movefile('meanafunctional_disc_4D.nii',procedure);
    movefile('wmeanafunctional_disc_4D.nii',procedure);
    movefile('swmeanafunctional_disc_4D.nii',procedure);
    movefile('rp_afunctional_disc_4D.txt',procedure);
end

end


