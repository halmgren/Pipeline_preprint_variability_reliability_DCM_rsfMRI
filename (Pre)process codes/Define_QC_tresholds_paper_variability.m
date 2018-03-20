function [threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_threshold_VOIs]=Define_QC_tresholds_paper_variability(procedure)

if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'ROI_Size')
    threshold_expl_var=60;
    threshold_max_conn=1/8;
    threshold_n_par_est=1;
    threshold_FD=1.5;
    threshold_threshold_VOIs=0.05;
end

end