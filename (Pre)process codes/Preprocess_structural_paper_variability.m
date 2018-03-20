function Preprocess_structural_paper_variability(dataset,number_subject,procedure,SPM_dir,Work_dir)

for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%02d',subject)]);
    session_number(1:2)=[];
    
    if strcmp(procedure,'Basic') %Same for all procedure, so only done once
        %only first session, because we use first anatomical image
        for session=1 %%leave it like this!!
            session_name=session_number(session).name;
            clear matlabbatch;
            
            disp(['subject: ' num2str(subject)]);
            disp('Preprocess_structural');
            
            cd(['sub-' sprintf('%02d',subject) '/' session_name '/anat']);
            [structural_images,~]=spm_select('FPList',pwd,'structural.nii');
            
            %Segmentation
            matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(structural_images);
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[SPM_dir '/tpm/TPM.nii,1']};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[SPM_dir '/tpm/TPM.nii,2']};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[SPM_dir '/tpm/TPM.nii,3']};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[SPM_dir '/tpm/TPM.nii,4']};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[SPM_dir '/tpm/TPM.nii,5']};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[SPM_dir '/tpm/TPM.nii,6']};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
            
            %Make bias-corrected skull-stripped image
            matlabbatch{2}.spm.util.imcalc.input(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
            matlabbatch{2}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
            matlabbatch{2}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
            matlabbatch{2}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
            matlabbatch{2}.spm.util.imcalc.output = 'Skullstr_biascor_structural';
            matlabbatch{2}.spm.util.imcalc.outdir = {''};
            matlabbatch{2}.spm.util.imcalc.expression = 'i1.*(i2+i3+i4)';
            matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{2}.spm.util.imcalc.options.mask = 0;
            matlabbatch{2}.spm.util.imcalc.options.interp = 1;
            matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
            
            %write structural to MNI space
            matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
            matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Skullstr_biascor_structural', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
            matlabbatch{3}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
            matlabbatch{3}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
            
            matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
            matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            spm_jobman('run',matlabbatch);
            
            mkdir(procedure);
            
            movefile('BiasField_structural.nii',procedure);
            movefile('c1structural.nii',procedure);
            movefile('c2structural.nii',procedure);
            movefile('c3structural.nii',procedure);
            movefile('c4structural.nii',procedure);
            movefile('c5structural.nii',procedure);
            movefile('mstructural.nii',procedure);
            movefile('Skullstr_biascor_structural.nii',procedure);
            movefile('structural_seg8.mat',procedure);
            movefile('wc2structural.nii',procedure);
            movefile('wc3structural.nii',procedure);
            movefile('wSkullstr_biascor_structural.nii',procedure);
            movefile('y_structural.nii',procedure);
            
        end
        
    end
    
end

end