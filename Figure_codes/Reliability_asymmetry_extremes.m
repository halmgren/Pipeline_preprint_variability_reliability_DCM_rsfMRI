function Reliability_asymmetry_extremes(SPM_dir,Work_dir)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This function uses an adapted version of scatterHistDiff
    %Downloaded from https://github.com/anne-urai/Tools/blob/master/plotting/scatterHistDiff.m
    %Adapted version is called scatterHistDiff_adapt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%
    %sub10
    %%%%%%%
    clear PP v c c_LEFT c_RIGHT;
    load([Work_dir '/DatasetKuehn/sub-08_summary/DCM/Basic/Smith/Full_model/GCM_DMN_full_estim.mat']);
    load([Work_dir '/DatasetKuehn/sub-08_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']);
    tmp=0;
    for session=1:length(GCM)

         if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
              disp(session);
              continue;
         end
         tmp=tmp+1;
         T = 0;
         C_LEFT = [0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0]'/3;
         C_RIGHT = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0]'/3;
         C= [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0]'/3;

         c_LEFT(tmp) = C_LEFT'*spm_vec(GCM{session}.Ep.A);
         c_RIGHT(tmp) = C_RIGHT'*spm_vec(GCM{session}.Ep.A);
         c(tmp) = C'*spm_vec(GCM{session}.Ep.A);
         v(tmp) = C'*GCM{session}.Cp(1:16,1:16)*C;
         PP(tmp)   = 1-spm_Ncdf(T,c(tmp),v(tmp));

    end
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    scatterHistDiff_adapt(c_LEFT,c_RIGHT,[],[],'k',PP);
    saveas(gcf,[Work_dir '/Figures_paper_variability/Reliability_asymmetry_sub10.bmp']);
    close;
    

    %%%%%%%
    %sub11
    %%%%%%%
    clear PP v c c_LEFT c_RIGHT;
    load([Work_dir '/DatasetGordon/sub-01_summary/DCM/Basic/Smith/Full_model/GCM_DMN_full_estim.mat']);
    load([Work_dir '/DatasetGordon/sub-01_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']);
    tmp=0;
    for session=1:length(GCM)

         if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
              disp(session);
              continue;
         end
        tmp=tmp+1;
         T = 0;
         C_LEFT = [0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0]'/3;
         C_RIGHT = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0]'/3;
         C= [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0]'/3;

         c_LEFT(tmp) = C_LEFT'*spm_vec(GCM{session}.Ep.A);
         c_RIGHT(tmp) = C_RIGHT'*spm_vec(GCM{session}.Ep.A);
         c(tmp) = C'*spm_vec(GCM{session}.Ep.A);
         v(tmp) = C'*GCM{session}.Cp(1:16,1:16)*C;
         PP(tmp)   = 1-spm_Ncdf(T,c(tmp),v(tmp));

    end

    figure('units','normalized','outerposition',[0 0 1 1]);
    scatterHistDiff_adapt(c_LEFT,c_RIGHT,[],[],'k',PP);
    saveas(gcf,[Work_dir '/Figures_paper_variability/Reliability_asymmetry_sub11.bmp']);
    close;
    clear;

end