function Plot_largest_connection(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%
%Biggest connection plotted against time
%Day2day dataset, subject 1
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function contains adapted version of spm_plot_ci_adapt (called spm_plot_ci_adapt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([Work_dir '/DatasetKuehn/sub-01_results/DCM/Basic/Smith/Full_model/PEB_A_mean_DMN.mat']);
load([Work_dir '/DatasetKuehn/sub-01_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']);
load([Work_dir '/DatasetKuehn/sub-01_summary/DCM/Basic/Smith/Full_model/GCM_DMN_full_estim.mat']);

flag=zeros(1,length(GCM));
for diagn=1:length(GCM)
    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
        flag(diagn)=1;
    end
end

GCM(find(flag==1))=[];

A_matrix=full(vec2mat(PEB.Ep,4,4))';
[a,b1]=max(A_matrix);
[a,b2]=max(max(A_matrix));
max_value=A_matrix(b1(b2),b2);

for session=1:length(GCM)
    Max_conn_strength(session)=GCM{session}.Ep.A(b1(b2),b2);
    Max_conn_var(session)=GCM{session}.Cp(b2*4+b1(b2),b2*4+b1(b2));
end

figure('units','normalized','outerposition',[0 0 1 1]);
spm_plot_ci_adapt(Max_conn_strength',full(diag(Max_conn_var)));
xlabel('Session','fontsize',36,'fontweight','bold');
ylabel('Posterior Expectation (Hz)','fontsize',36,'fontweight','bold');
title('Connection from lIPC to PRC','fontsize',36,'fontweight','bold')
set(gca,'fontsize',30,'fontweight','bold');
set(gcf,'PaperPositionMode','auto');

saveas(gcf,[Work_dir '/Figures_paper_variability/Strongest_connection_S3.bmp']);
close;

end