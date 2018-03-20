function BS_violin_plots(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%
%This code uses an adapted version of rst_data_plot
%downloaded from https://github.com/CPernet/Robust_Statistical_Toolbox/blob/master/graphic_functions/rst_data_plot.m
%adapted version is called rst_data_plot_adapt
%%%%%%%%%%%%%%%%%%
ntwrk_name='DMN';
procedure='Basic';
name_ROI_def='Smith';

load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BS_correlation_' ntwrk_name '.mat'],'BS_correlation');
load([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BS_correlation_delat_' ntwrk_name '.mat'],'BS_correlation_delat');

rst_data_plot_adapt([BS_correlation; BS_correlation_delat]','estimator','mean','point_size',80);
hold on;

set(gca,'fontsize',30,'fontweight','bold');
set(gca,'XTick',[1 2.25],'XTickLabel',{'Original','Dominance-ordered'});

%label
yL=ylabel('Correlation','Fontweight','bold');

%axes
ax=gca;
yrule=ax.YAxis;
yrule.FontSize=30;
xrule=ax.XAxis;
xrule.FontSize=38;
yL.FontSize=50;


set(gcf,'PaperPositionMode','auto');        
saveas(gcf,[Work_dir '/Figures_paper_variability/Figure_BS_violin_plots.bmp']);
close;

end